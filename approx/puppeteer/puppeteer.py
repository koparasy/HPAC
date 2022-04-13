import argparse
import shutil
import time
import json
import sys
import os
import yaml
from pathlib import Path
from yaml import CLoader
import subprocess
import glob
import math
from string import Template
import re
from SALib.sample import saltelli as sobol_sample
from SALib.sample import morris as morris_sample
from SALib.analyze import sobol as sobol_analyze
from SALib.analyze import morris as morris_analyze
import pandas as pd
import numpy as np
from enum import IntEnum
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

tex_fonts = {
# Use LaTeX to write all text
"text.usetex": True,
"font.family": "serif",
# Use 10pt font in plots, to match 10pt font in document
"axes.labelsize": 8,
"font.size": 8,
# Make the legend/label fonts a little smaller
"legend.fontsize": 6,
"xtick.labelsize": 8,
"ytick.labelsize": 8
}

def set_size(width, fraction=1, sH = 1):
    """Set figure dimensions to avoid scaling in LaTeX.
    Parameters
    ----------
    width: float
    Document textwidth or columnwidth in pts
    fraction: float, optional
    Fraction of the width which you wish the figure to occupy
    Returns
    -------
    fig_dim: tuple
    Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27
    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2
    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = sH *fig_width_in * golden_ratio
    fig_dim = (fig_width_in, fig_height_in)
    return fig_dim

plt.rcParams.update(tex_fonts)

textWidth=505.89


def closest_power(x):
  possible_results = math.floor(math.log(x, 2)), math.ceil(math.log(x, 2))
  return 2**(min(possible_results, key= lambda z: abs(x-2**z)))

class PetrubationMethod(IntEnum):
  ABS=1
  REL=2
  DIST=3

def revertHPACToSaLib(problem, index_to_keys, input_space_dict, samples):
  Y=[]
  X = []
  for i,vals in enumerate(samples):
    errors = list()
    for j,(label,varName) in enumerate(index_to_keys):
      if isinstance(problem, np.ndarray):
        errors.append(problem[i][j])
      else:
        errors.append(vals['exp'][label][varName][1])
    X.append(errors)
    Y.append(float( vals['Error'][0]))
  return X,Y

def vars_sample(index_to_keys, vars_input_file, exponents=True):
  extract_data = False
  error_vals = []
  with open(vars_input_file) as sampled_file:
    for line in sampled_file:
      if extract_data:
        vals = [float(s) for s in line.split()]
        if len(index_to_keys) == len(vals):
          error_vals.append(vals)
        else:
          raise ValueError("Lengths of keys and values are different", len(index_to_keys), len(vals))
      else:
        if 'factor' in line:
          extract_data = True
  return np.array(error_vals)

def convertProblemToHPAC(problem, index_to_keys, input_space_dict, error_method):
  record = []
  for i,vals in enumerate(problem):
    tmp = dict()
    for k,v in enumerate(vals):
      label, varName = index_to_keys[k]
      if label not in tmp:
        tmp[label] = dict()
      if varName not in tmp[label]:
        tmp[label][varName] = input_space_dict[label][varName].copy()
      tmp[label][varName][1] = v
      tmp[label][varName][3] = int(error_method)
    experiment = {}
    experiment['exp'] = tmp
    experiment['Status'] = 'Not Run'
    record.append(experiment)
  return record

def createProblemDescr(space, error_method, error_descr):
  names = []
  bounds = []
  index_to_keys = []
  for label, outputs in space.items():
    if label in error_descr:
      lb = error_descr[label]['lower_bound']
      ub = error_descr[label]['upper_bound']
    else:
      lb = error_descr['all']['lower_bound']
      ub = error_descr['all']['upper_bound']
    for varName, values in outputs.items():
      if len(outputs.items()) != 1:
        names.append(f'{label}_{varName}')
      else:
        names.append(f'{label}')
      bounds.append([lb, ub])
      index_to_keys.append((label, varName))
  problem = {
    'num_vars' :  len(names),
    'names':  names,
    'bounds' : bounds
  }
  return (index_to_keys, problem)

def createWorkerCMD(args, chunk):
  cmd = f'python {sys.argv[0]}'
  cmd += ' -i %s' % args.input
  cmd += ' -r %s' % args.root
  cmd += ' -t %s' % args.sample_technique
  cmd += ' -m %s' % args.method
  cmd += ' -c %s' % args.cluster
  cmd += ' -s %d' % chunk
  cmd += ' -e run'
  cmd += ' -w $worker_id'
  return cmd

def TransformExponentToValueSpace(expSpace, sample_technique, index_to_keys, error_descr):
  for index,value in enumerate(index_to_keys):
    if value[0] in error_descr:
      key=value[0]
    else:
      key = 'all'
    print(f'Key is {key} for {index}')
    error_direction=error_descr[key]['error_direction']
    base =error_descr[key]['base']
    max_error = error_descr[key]['max_error']
    upper_bound = error_descr[key]['upper_bound']
    exponent_space = error_descr[key]['exponent_space']

    if exponent_space:
      if error_direction == 'Positive':
        expSpace[:,index] =np.abs(expSpace[:,index]) * np.full(expSpace[:,index].shape,base)**(np.abs(expSpace[:,index])*(-1.0) + max_error)
      elif error_direction == 'Negative':
        expSpace[:,index]= (-1.0 * np.abs(expSpace[:,index])) * np.full(expSpace[:,index].shape,base)**(np.abs(expSpace[:,index])*(-1.0)+max_error)
      elif error_direction == 'Both':
        if sample_technique == 'vars':
          expSpace[:,index] = np.sign(expSpace[:,index]) * np.full(expSpace[:,index].shape ,base)**(np.abs(expSpace[:,index])*(1.0) + max_error - upper_bound)
        else:
          expSpace[:,index] = np.sign(expSpace[:,index]) * np.full(expSpace[:,index].shape ,base)**(np.abs(expSpace[:,index])*(-1.0) + max_error)
    else:
      expSpace[:,index] = expSpace[:,index]
  return expSpace

def submitReduction(log_dir, jobDir, dependencies):
    tmp =       ('#!/bin/bash',
                '#SBATCH --time 01:59:00',
                '#SBATCH --open-mode append',
                f'#SBATCH --output {log_dir}/reduction.out',
                f'#SBATCH --error {log_dir}/reduction.err',
                'python $CMD'
                )
    cmd = ' '.join(sys.argv)
    sbatchTMP = '\n'.join(tmp)
    template = Template(sbatchTMP)
    script = template.substitute(CMD=cmd)
    batch_file = f'{jobDir}/master.batch'
    with open(batch_file,'w') as fd:
        fd.write(script)
    return submitJob(batch_file, 1, 1, dependencies)

def createBatchScript(node_dir, node_id, workers_dirs, cmd, time='00:30:00'):
  worker_dirs_str=' '.join(workers_dirs)
  with open('batch_templates/distributeWork.batch', 'r') as fd:
    template = Template(fd.read())
  script = template.substitute(TIME=time, NODE_ID=node_id, OUTDIR=node_dir, WORKER_IDS=worker_dirs_str, CMD=cmd)
  with open(f'{node_dir}/{node_id}.batch', 'w') as fd:
    fd.write(script)
  return f'{node_dir}/{node_id}.batch'

def submitJob(script, num_cores, cores_per_task, dependencies=None):
    submitCmd = f'sbatch -N 1 -n {num_cores} -c {cores_per_task}'
    if dependencies != None:
        submitCmd += f' {dependencies}'
    submitCmd +=  f' {script}'
    print(submitCmd)
    p = subprocess.run( submitCmd, shell=True, capture_output=True)
    out=p.stdout.decode('utf-8')
    err=p.stderr.decode('utf-8')
    vals =  re.compile('Submitted batch job (\d+)').findall(out)
    if len(vals) != 1:
        print('Could Not Match Everything')
        print(out)
        print(err)
        sys.exit()
    jobId = int(vals[0])
    return jobId

def computeNumJobs(numExperiments, cluster):
  print(cluster)
  NumTasks = int(math.floor(cluster['NumCores']/cluster['CoresPerTask']))
  TotalWorkers = NumTasks * cluster['NumNodes']
  ExperimentsPerTask = math.ceil(numExperiments / TotalWorkers)
  return cluster['NumNodes'], TotalWorkers, ExperimentsPerTask, NumTasks

def updateResults(experiments, results):
  analyzed = 0
  for index, val in enumerate(experiments):
    if 'Error' in val:
      analyzed+=1
    else:
      fn=f'{results}/exp_{index}.json'
      if Path(fn).is_file():
        with open(fn,'r') as fd:
          result = json.load(fd)
        experiments[index]['Status'] = 'Finished'
        experiments[index]['Runtime'] = result['Runtime']
        experiments[index]['Error'] = result['Error']
        analyzed += 1
    if analyzed % 1000 == 1:
      print('Analyzed ', analyzed, 'Out of ', len(experiments))
  return analyzed, experiments

def createDir(path):
    testPath = Path(path)
    if testPath.exists():
        if not testPath.is_dir():
            print(f'Path already exists {path} and it is not a directory')
            sys.exit()
    else:
        testPath.mkdir(parents=True)
    return path

def compile(bench, config):
  build_dir = config['build_dir']
  for k, v in config['versions'].items():
    p = subprocess.run( v['clean'], cwd=build_dir, shell=True)
    if p.returncode != 0:
      print('could not clean up {build_dir}--{k}')
      sys.exit()
    p = subprocess.run( v['build'], cwd=build_dir, shell=True)
    if p.returncode != 0:
      print('could not clean up {build_dir}--{k}')
      sys.exit()

def computeQuality(config, accurate_path, test_path, stdout=None, stderr=None):
  cmd = '%s/build/bin/quality -m %s -a %s -t %s' % (config['root_dir'], config['metric'], accurate_path, test_path)
  p = subprocess.run( cmd, shell=True, capture_output=True)
  out = p.stdout.decode('utf-8')
  err = p.stderr.decode('utf-8')
  if stdout != None:
    with open (stdout, 'a') as fd:
      fd.write(out)

  if stderr != None:
    with open (stderr, 'a') as fd:
      fd.write(err)
  return out,err

def execute(bench, config, version, outputName, stdout=None, stderr=None, env=None):
  cmd = '%s %s' % (config['versions'][version]['exe'], config['cmd'])
  cmd = cmd.replace('<input>', config['input'][0])
  cmd = cmd.replace('<output>', outputName)
  print(cmd)
  p = subprocess.run( cmd, cwd=config['build_dir'], shell=True, capture_output=True, env=env)
  out = p.stdout.decode('utf-8')
  err = p.stderr.decode('utf-8')
  if stdout != None:
    with open (stdout, 'a') as fd:
      fd.write(out)

  if stderr != None:
    with open (stderr, 'a') as fd:
      fd.write(err)
  return out,err

def main():
    parser = argparse.ArgumentParser(description='Create Sample Space for Sobol')
    parser.add_argument('-i', '--input', dest='input', type=str, help='Yaml file containing the description of the application', required=True)
    parser.add_argument('-r', '--root-dir', dest='root', type=str, help='Output directory containing all outpus produced by experiments', required=True)
    parser.add_argument('-t', '--technique', dest='sample_technique', type=str.lower, choices=('vars', 'sobol', 'morris'), help='Technique to use to sample space', required=True)
    parser.add_argument('-e', '--execution-type', dest='type', choices=('setup', 'deploy', 'run', 'analyze', 'visualize', 'clean'), help='Execution type', required=True)
    parser.add_argument('-c', '--cluster-descr', dest='cluster', help='cluster description', required=True)
    parser.add_argument('-w', '--worker-id', dest='wid', help='Id of the current worker')
    parser.add_argument('-s', '--size-work', dest='sw', help='Items required to be checked by worker')
    parser.add_argument('-v', '--vars-input', dest='vars', type=str, help='input file containing the samples of the VARS technique')
    parser.add_argument('-m', '--petrubation-method', type=str.upper, dest='method', choices=('ABS', 'REL', 'DIST'), help='Method to inject errors', required=True)
    parser.add_argument('-p', '--petrubation-description', type=str, dest='petr_descr', help='description of petrubation')
    parser.add_argument('-N', '--number-of-samples', type=int, dest='num_samples', help='Number Of Samples')

    args = parser.parse_args()

    rootDir = os.path.abspath(createDir(args.root))
    dbDir = createDir(f'{rootDir}/{args.method}/{args.sample_technique}')
    interm = createDir(f'{dbDir}/interm')
    logDir = createDir(f'{dbDir}/logs')
    nodeDir = createDir(f'{dbDir}/nodes')
    expDir = createDir(f'{dbDir}/experiments/')
    results = createDir(f'{dbDir}/results')
    sampleDir = createDir(f'{dbDir}/sample')
    error_method = PetrubationMethod[args.method]

    sample_file = f'{dbDir}/sample/space.json'
    space_file=f'{dbDir}/space.json'

    with open(args.input, 'r') as f:
        config = yaml.load(f, Loader=CLoader)
    #This is a hack. 
    bench = list(config.keys())[0]
# Setup: Sets an experimentation campaign. Creates directories and samples for our function envaluation. Pe
    if args.type == 'setup':
      if args.petr_descr is None:
        print('Setup requires a petrubation yaml file to describe how to transform sampled space to petrubation space')
        sys.exit()

      if args.num_samples is None:
        print('Setup requires to define the number of samples')
        sys.exit()

      if (args.sample_technique == 'vars') and (args.vars is None):
        print('The vars technique requires option -v/--vars-input')
        sys.exit()

      if args.sample_technique == 'vars':
        shutil.copyfile(args.vars, f'{dbDir}/vars_orig_input.txt')

      # I need to compile everything to make sure 
      # binaries exist before deploying
      compile(bench, config[bench])
      # I need to create accurate output file
      execute(bench, config[bench], 'accurate', f'{interm}/accurate.out', f'{logDir}/accurate.out', f'{logDir}/accurate.err')
      current_env = os.environ.copy()
      current_env['OMP_NUM_THREADS'] = str(1)
      current_env['PETRUBATE_TYPE'] = 'RECORD'
      current_env['PETRUBATE_FILE'] = space_file
      execute(bench, config[bench], 'petrubate', f'{interm}/record.out', f'{logDir}/record.out', f'{logDir}/record.err', current_env)
      with open(space_file, 'r') as fd:
        input_space_dict = json.load(fd)

      with open(args.petr_descr, 'r') as fd:
        error_descr = yaml.load(fd, Loader=CLoader)

      index_to_keys, problemDescr = createProblemDescr(input_space_dict, error_method, error_descr)
      print(index_to_keys)

      with open(f'{dbDir}/index_to_keys.json', 'w') as fd:
        json.dump(index_to_keys, fd)
      with open(f'{dbDir}/problem_descr.json', 'w') as fd:
        json.dump(problemDescr, fd)

      if args.sample_technique == 'sobol':
        number_parameter = problemDescr['num_vars']
        N = args.num_samples
        print('Before finding power of 2', N)
        N = closest_power(N)
        print('After finding power of 2', N)
        param_values = sobol_sample.sample(problemDescr, N, calc_second_order=True)
      elif args.sample_technique == 'morris':
        number_parameter = problemDescr['num_vars']
        N= args.num_samples
        optimal_trajectories=70
        num_levels = int(N)
        num_levels = num_levels - num_levels % 2
        print('N IS:', N)
        param_values = morris_sample.sample(problemDescr, N=N, seed=0, local_optimization = True, optimal_trajectories=optimal_trajectories)
#        param_values = morris_sample.sample(problemDescr, N=N, seed=0, local_optimization = True, num_levels=num_levels)
      elif args.sample_technique == 'vars':
        param_values =  vars_sample(index_to_keys, f'{dbDir}/vars_orig_input.txt')
      print(param_values.shape)

      np.save(f'{dbDir}/sampled_space.npy', param_values)

      param_values = TransformExponentToValueSpace(param_values,args.sample_technique,index_to_keys, error_descr)

      for i, val in enumerate(index_to_keys):
        tmp = param_values[:,i]
        print(val[0], 'Threshold for Min Negative val',np.min(tmp[tmp<0], axis=0))
        print(val[0], 'Threshold for Max Negative val',np.max(tmp[tmp<0], axis=0))
        print(val[0], 'Threshold for Min Positive val',np.min(tmp[tmp>0], axis=0))
        print(val[0], 'Threshold for Max Positive val',np.max(tmp[tmp>0], axis=0))

      with open(f'{dbDir}/error_descr.yaml', 'w') as fd:
        yaml.dump(error_descr,fd, Dumper=yaml.CDumper, default_flow_style=False)

      np.save(f'{dbDir}/transformed_sampled_space.npy', param_values)

      print(param_values.max(),param_values.min())
      experiments = convertProblemToHPAC(param_values, index_to_keys, input_space_dict, error_method)
      print(param_values.shape)
      with open(sample_file,'w') as fd:
        json.dump(experiments,fd)
      print(f'Created {len(experiments)} experiments, Problem:{len(param_values)}')
    elif args.type == 'deploy':
      with open(sample_file, 'r') as fd:
        experiments = json.load(fd)
      analyzed, experiments = updateResults(experiments, results)
      with open(sample_file, 'w') as fd:
        json.dump(experiments,fd)

      print(f'Number of total Experiments: {len(experiments)} Currently Analyzed: {analyzed}')
      if len(experiments) == analyzed:
        return

      with open(args.cluster, 'r') as fd:
        cluster = yaml.load(fd, Loader=CLoader)
      numJobs, TotalWorkers, ExperimentsPerTask, WorkersPerJob = computeNumJobs(len(experiments), cluster)
      print( numJobs, TotalWorkers, ExperimentsPerTask )
      jobs = []
      cmd = createWorkerCMD(args, ExperimentsPerTask)
      for i in range(numJobs):
        worker_ids = [str(k) for k in range(i*WorkersPerJob, (i+1)*WorkersPerJob)]
        script = createBatchScript(nodeDir, f'node_{i}', worker_ids, cmd)
        job = submitJob(script, WorkersPerJob, cluster['CoresPerTask'])
        jobs.append(job)
      deps=','.join([str(j) for j in jobs])
      dependencies = f'--dependency=afterany:{deps}'
      submitReduction(nodeDir,nodeDir, dependencies)
    elif args.type == 'run':
      wid = int(args.wid)
      wSize = int(args.sw)

      with open(sample_file, 'r') as fd:
        experiments = json.load(fd)

      with open(args.cluster, 'r') as fd:
        cluster = yaml.load(fd, Loader=CLoader)

      tempDir = cluster['localDirectory']

      start = wid * wSize
      end = min((wid+1)*wSize, len(experiments))
      print(start,end)
      current_env = os.environ.copy()
      current_env['OMP_NUM_THREADS'] = str(1)
      current_env['PETRUBATE_TYPE'] = 'PETRUBATE'
      for i in range(start,end):
        start = time.time()
        exp = experiments[i]
        if exp['Status'] != 'Not Run':
          continue
        errorDescr = exp['exp']
        expName = f'{tempDir}/descr_{wid}.json'
        with open(expName, 'w') as fd:
          json.dump(errorDescr, fd)
        current_env['PETRUBATE_FILE'] = expName
        out, err = execute(bench, config[bench], 'petrubate', f'{tempDir}/qoi_{wid}.out', f'{logDir}/exp_{wid}.out', f'{logDir}/exp_{wid}.err', current_env)
        items = re.search(config[bench]['measure'],out)
        if items != None:
          rt = list(items.groups())
          runtime = [t for t in rt if t]
        else:
          runtime=-1
        out, err = computeQuality(config[bench], f'{interm}/accurate.out', f'{tempDir}/qoi_{wid}.out', f'{logDir}/quality_{wid}.out', f'{logDir}/quality_{wid}.err')
        items = re.search(config[bench]['quality_pattern'],out)
        qual= list(items.groups())
        quality = [q for q in qual if q]
        with open(f'{results}/exp_{i}.json', 'w') as fd:
          json.dump({'Runtime' : runtime , 'Error' : quality }, fd)
        end = time.time()
        print(f'Iter:{i} Took {end - start} seconds')
    elif args.type == 'clean':
      with open(sample_file, 'r') as fd:
        experiments = json.load(fd)
      numExperiments = len(experiments)
      for i in range (numExperiments):
        expName = f'{results}/exp_{i}.json'
        try:
          os.remove(expName)
        except:
          print("Error while deleting file ", expName)
        print(f'Deleted file {expName}')

    elif args.type == 'analyze':

      with open(space_file, 'r') as fd:
        input_space_dict = json.load(fd)

      with open(sample_file, 'r') as fd:
        samples = json.load(fd)

      with open(f'{dbDir}/index_to_keys.json', 'r') as fd:
        index_to_keys = json.load(fd)

      with open(f'{dbDir}/problem_descr.json', 'r') as fd:
        problemDescr = json.load(fd)

      param_values = np.load(f'{dbDir}/sampled_space.npy')
      X,Y = revertHPACToSaLib(param_values, index_to_keys, input_space_dict, samples)
      Y = np.array(Y)
      X = np.array(X)

      np.save(f'{dbDir}/application_evaluation.npy', Y)

      print('X-Shape:', X.shape, 'Y Shape:', Y.shape)
      if args.sample_technique == 'vars':
        result = np.hstack((X, np.atleast_2d(Y).T))
        np.savetxt(f'{dbDir}/vars_out.txt', result)
        with open(f'{dbDir}/names.txt', 'w') as fd:
          for i,v in enumerate(index_to_keys):
            fd.write(f'{i}=>{v}\n')
        sys.exit()

      sensitivities = {}
      num_resamples = 1000
      if args.sample_technique == 'sobol':
        Si = sobol_analyze.analyze(problemDescr, Y, print_to_console=True, num_resamples=num_resamples)
        sensitivities['Sensitivity'] = list(Si['ST'])
        sensitivities['Confidence'] = list(Si['ST_conf'])
      elif args.sample_technique == 'morris':
        num_levels = int(0.1*len(Y))
        print(num_levels)
        Si = morris_analyze.analyze(problemDescr, X, Y, conf_level=0.95, print_to_console=True, num_levels=num_levels, seed=0, num_resamples=num_resamples)
        sensitivities['Sensitivity'] = list(Si['mu_star'].filled())
        sensitivities['Confidence'] = list(Si['mu_star_conf'])
      sensitivities['names']  = problemDescr['names']

      with open(f'{dbDir}/analysis.json', 'w') as fd:
        json.dump(sensitivities, fd)
    elif args.type == 'visualize':
      viz_dir = createDir(f'{dbDir}/Figures/')
      with open(f'{dbDir}/index_to_keys.json', 'r') as fd:
        index_to_keys = json.load(fd)

      regions = [ name[0].replace('_','') for name in index_to_keys ]
      for v in ['sampled_space', 'transformed_sampled_space']:
        X = np.load(f'{dbDir}/{v}.npy')
        columns = X.shape[1]
        for i in range(columns):
          x = X[:,i]
          print(x.shape)
          df = pd.DataFrame(x,columns=[index_to_keys[i][0]])
          sizes=set_size(width=textWidth, fraction=0.5)
          fig, ax = plt.subplots(figsize=sizes)
          ax = df.plot.hist(bins=100, alpha=0.5, ax = ax)
          ax.figure.savefig(f'{viz_dir}/{v}_{index_to_keys[i][0]}_dist.pdf',bbox_inches='tight')
          plt.close()
      Y = np.load(f'{dbDir}/application_evaluation.npy')
      df = pd.DataFrame(Y,columns=['Application Output'])
      sizes=set_size(width=textWidth, fraction=0.5)
      fig, ax = plt.subplots(figsize=sizes)
      ax = df.plot.hist(bins=100, alpha=0.5, ax = ax)
      ax.figure.savefig(f'{viz_dir}/application_out.pdf',bbox_inches='tight')
      plt.close()

main()

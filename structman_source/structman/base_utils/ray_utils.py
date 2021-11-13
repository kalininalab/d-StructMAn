import ray
import os
from structman import settings
import subprocess

def ray_init(config):
    if ray.is_initialized():
        return

    os.environ["PYTHONPATH"] = f'{settings.ROOT_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.LIB_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.RINERATOR_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.OUTPUT_DIR}:{os.environ.get("PYTHONPATH", "")}'
    if config.iupred_path != '':
        os.environ["PYTHONPATH"] = f'{os.path.abspath(os.path.realpath(config.iupred_path))}:{os.environ.get("PYTHONPATH", "")}'

    logging_level = 20
    if config.verbosity <= 1:
        logging_level = 0
    ray.init(num_cpus=config.proc_n, include_dashboard=False, ignore_reinit_error=True, logging_level = logging_level)

def ray_hack():
    # hack proposed by the devs of ray to prevent too many processes being spawned
    return
    #resources = ray.ray.get_resource_ids()
    #cpus = [v[0] for v in resources['CPU']]
    #psutil.Process().cpu_affinity(cpus)

def monitor_ray_store():
    p = subprocess.Popen(['ray','memory'])

    page, err = p.communicate()
    print(page)

    return

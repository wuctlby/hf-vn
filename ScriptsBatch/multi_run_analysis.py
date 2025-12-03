import argparse
import yaml
import copy
import os
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
RUN_ANALYSIS_SCRIPT = os.path.join(SCRIPT_DIR, '..', 'run_analysis.py')

meta_task = []

def compose_mode(task_config, mode):
    operations = task_config["operations"]
    for operation in operations.keys():
        if operation in mode["operations"]:
            operations[operation] = mode["operations"][operation]
        else:
            operations[operation] = False
    return task_config

def compose_task(itask, task, mode):
    global meta_task
    with open(task["config"], 'r') as f:
        config = yaml.safe_load(f)
    task_config = copy.deepcopy(config)
    task_config = compose_mode(task_config, mode)
    for key, value in task.items():
        if key not in meta_task:
            if key == "dataDir":
                task_config["preprocess"]["data"]["2023"]["files"] = value
            elif key == "resoDir":
                task_config["preprocess"]["data"]["2023"]["resolution"] = value
            else:
                task_config[key] = value

    composed_task_path = f"{os.path.dirname(task['config'])}/composed_{itask}_{os.path.basename(task['config'])}"
    with open(composed_task_path, 'w') as f:
        yaml.dump(task_config, f)

    if task['method'] == 'correlated':
        cmd = (
            f"python {RUN_ANALYSIS_SCRIPT} "
            f"{composed_task_path} -corr"
        )
    else:
        cmd = (
            f"python {RUN_ANALYSIS_SCRIPT} "
            f"{composed_task_path} -comb"
        )
    os.system(cmd)
    os.remove(composed_task_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config_file", "-c", type=str,
                        default="config_multi_run_analysis.yml", help="input config file")
    args = parser.parse_args()

    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)

    tasks = [config[task] for task in config if task.startswith("task")]
    modes = [config[mode] for mode in [task["mode"] for task in tasks]]
    meta_task.extend([meta for meta in config["META_TASK"].keys()])
    
    for itask, (task, mode) in enumerate(zip(tasks, modes)):
        compose_task(itask, task, mode)
import os

def mkdir_if_not_exists(default_save_path: str):
    if not os.path.exists(default_save_path):
        os.mkdir(default_save_path)


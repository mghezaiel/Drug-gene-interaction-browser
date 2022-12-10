import os 
from app import App

# DIRS
DIR = os.path.dirname(os.getcwd())
DATA_DIR = os.path.join(DIR,"data")
CACHE_DIR = os.path.join(DIR,"cache")
    
if __name__ == "__main__":
    
    # LAUNCH THE APP
    App(data_dir = DATA_DIR, 
        cache_dir = CACHE_DIR)
    
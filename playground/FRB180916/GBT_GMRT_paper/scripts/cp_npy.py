import glob
import os
import shutil
import sys

src = sys.argv[1]
dest = sys.argv[2]

def copy():
    for file_path in glob.glob(os.path.join(src, '**', '*.npy'), recursive=True):
        new_path = os.path.join(dest, os.path.basename(file_path))
       	shutil.copy(file_path, new_path)
copy()

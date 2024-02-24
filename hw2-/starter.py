import os

os.environ["LD_LIBRARY_PATH"]="/act/opt/intel/composer_xe_2013.3.163/mkl/lib/intel64:"+os.environ.get("LD_LIBRARY_PATH","")
os.system("make")
os.system("./main lapack 4096")
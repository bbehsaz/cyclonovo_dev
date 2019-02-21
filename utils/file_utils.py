import os
import errno
import sys
from os.path import isdir, isfile, getsize


def error(msg):
    print(msg)
    sys.exit(1)


def verify_file(fpath, description=''):
    if not check_file(fpath):
        error(fpath + (' (' + description + ' file)' if description else '') + ' is empty or does not exist!')


def verify_dir(dirpath, description=''):
    if not isdir(dirpath):
        error(dirpath + (' (' + description + ' dir)' if description else '') + ' does not exist!')


def check_file(fpath):
    return isfile(fpath) and getsize(fpath) > 0


def remove_if_exists(fpath):
    # http://stackoverflow.com/questions/10840533/most-pythonic-way-to-delete-a-file-which-may-not-exist
    try:
        os.remove(fpath)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


def count_lines(fpath):
    with open(fpath) as f:
        return len(f.readlines())


def removeDir(path):
    #https://stackoverflow.com/questions/303200/how-do-i-remove-delete-a-folder-that-is-not-empty-with-python
    import os
    import stat
    import shutil
    def errorRemoveReadonly(func, path, exc):
        excvalue = exc[1]
        if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
            # change the file to be readable,writable,executable: 0777
            os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  
            # retry
            func(path)
        # else:
        #     # raiseenter code here
    shutil.rmtree(path, ignore_errors=False, onerror=errorRemoveReadonly) 



def mkdir_p(path):
    import errno
    import os
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
import xml.etree.ElementTree
from distutils import dir_util
from os.path import join, isfile, isdir, basename


def parse_params_xml(fpath):
    with open(fpath) as f:
        content = f.read()
    params = dict()
    file_mapping = dict()
    for e in xml.etree.ElementTree.fromstring(content).findall('parameter'):
        if e.attrib['name'] == 'upload_file_mapping':
            fname, real_fname = e.text.split('|')[0:2]
            file_mapping[basename(fname)] = real_fname
            continue
        if e.text == 'on':
            value = True
        elif e.text == 'off':
            value = False
        else:
            value = e.text
        params[e.attrib['name']] = value
    return params, file_mapping



def removeDir(path):
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
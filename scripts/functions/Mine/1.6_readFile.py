def readFile(path):
    '''
    reads file to var
    Arguments:
        file path
    Vars- 
        f: opened file
        file: read file
    Return:
        file
    '''
    f=open(path)
    file= f.read()
    f.close()
    return file
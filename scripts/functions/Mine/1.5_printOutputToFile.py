import sys

def printOutputToFile(path, type, output):
    '''
    temporarily redirects stdout to a file.
    Arguments:
        takes a path to output the file and its type (r, w, a) for example. And the output to print.
    Vars -
        stdout: saves unchanged standard output value
    Returns:
        None
    '''
    stdout= sys.stdout
    sys.stdout= open(path, type)
    print(output)
    sys.stdout= stdout
    return None




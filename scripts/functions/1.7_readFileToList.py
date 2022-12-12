
def readFileToList(path_to_file: str):
    '''
    reads a file line by line. Adds lines to a list and removes repeated and empty items.
    Args:
        takes the path to the file as an argument
    Vars - 
        fileList: unchanged file content
        List: transformation of the file content into a python list.
    Returns:
        cleared and transformed list.
    '''
    with open(path_to_file, 'r') as fileList:
        List = [List.rstrip() for List in fileList]
    List= list(set(List)) 
    if "" in List:
        List.remove("")
    return List



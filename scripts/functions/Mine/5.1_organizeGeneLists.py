import os


def organizeSpecieGeneList(specie):
    '''
    Organizes the gene list and reacts incase list is empty.
    Arguments:
        takes the specie name as an argument
    Vars-
        None
    Returns:
        None
    '''
    if os.path.exists('./GeneLists'):
        os.system('rm ./GeneLists/*')
    else: os.system('mkdir ./GeneLists')

    try:
        os.rename('GeneList.txt', './GeneLists/{}_GeneList'.format(specie))
    except:
        f= open('./GeneLists/noGeneSpecies.txt', 'a')
        f.write(specie +'\n')
        print(f'The following species has zero genes: {specie}')
    return None
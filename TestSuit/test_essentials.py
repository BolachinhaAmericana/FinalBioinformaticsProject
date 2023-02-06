import sys
from scripts.essentials import get_user_arguments, entrez_search, entrez_fetch, print_to_file, file_reader
from Bio import Entrez

def test_get_user_arguments():
    number_of_args = 2
    saved_argv = sys.argv
    sys.argv = ['script_name', 'arg1', 'arg2']
    result = get_user_arguments(number_of_args)
    assert result == ['arg1', 'arg2']
    sys.argv = saved_argv

def test_entrez_search():
    database = 'Nucleotide'
    term = 'test'
    result = entrez_search(database, term)
    assert isinstance(result, dict)
    assert 'IdList' in result
    
    database = 'Gene'
    result = entrez_search(database, term)
    assert isinstance(result, dict)
    assert 'IdList' in result
    
    database = 'Protein'
    result = entrez_search(database, term)
    assert isinstance(result, dict)
    assert 'IdList' in result

def entrez_fetch(database, web_environment, query_key):
    try:
        fetch_handle = Entrez.efetch(db=database,
                                     retmode='text',
                                     webenv=web_environment,
                                     query_key=query_key)
        result = fetch_handle.read()
        fetch_handle.close()
        return result
    except HTTPError as e:
        print(f'HTTP error: {e}')
        return None

def test_print_to_file():
    output_path = 'test_output.txt'
    mode = 'w'
    output = 'test'
    result = print_to_file(output_path, mode, output)
    assert result == None
    
    with open(output_path, 'r') as file:
        content = file.read()
        assert content == 'test\n'
        
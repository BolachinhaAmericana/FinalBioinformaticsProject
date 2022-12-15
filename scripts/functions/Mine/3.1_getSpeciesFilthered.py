
# Not a function!!

path = 1

f=open(path)
file= f.read()
f.close()

filther= file.replace('<ScientificName>', '\n')
# Send output to filthering.txt
# subprocess.call('./bashScripts/getSpeciesFilthered.sh')
# output = loading.txt

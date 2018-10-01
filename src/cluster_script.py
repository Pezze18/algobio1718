#!/usr/bin/env python3
import os
import time
import datetime
import paramiko
import getpass
import lib.inout as inout
from lib.inout import *

"""
    --- WARNING ---
    FIRST TIME INFO, MANDATORY TODOs

    (1)     Se non hai mai avviato questo script, come prima cosa esegui il login:
            >>> ssh username@login.dei.unipd.it
            >>> git clone <url>
            clonare la repository dentro la 'home' del nostro spazio sul DEI.

    (2)     Per evitare inutili disagi, eseguire questo script sempre da dentro la cartella src/.
"""

#ho importato i parametri di default da inout.py
parameters['proteins_input']="../"+parameters['proteins_input']
parameters['samples_input']="../"+parameters['samples_input']
parameters['genes_input']="../"+parameters['genes_input']
parameters['filter_input']="../"+parameters['filter_input']

server = "login.dei.unipd.it"

# Crea un'istanza del client SSH
ssh = paramiko.SSHClient()
# <<Required, since "login.dei.unipd.it" is not a "well known ssh host">> (???)
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Per evitare di scrivere il proprio username, se la macchina dalla quale lo sia avvia è "nota"
if getpass.getuser() == 'DavideDP':
    pwd = str(getpass.getpass(prompt="Password: "))
    username = "dallepez"
elif getpass.getuser() == 'iltuonomesultuoPC':
    pwd = str(getpass.getpass(prompt="Password: "))
    username = "iltuonomeDEI"
else:
    username = str(input("Please, enter your username (@" + server + "): "))
    pwd = str(getpass.getpass(prompt=username + "@" + server + ": "))

# Connesione al server
ssh.connect(server, username=username, password=pwd)

# Apertura di un'istanza per il file transfer (con secure file transfer protocol)
sftp = ssh.open_sftp()  # per trasferire i file

# Cartella remota del nostro progetto
remote_path = "/home/" + username + "/algobio1718/"

# Cartella locale del nostro progetto
local_path = os.path.dirname(os.getcwd()) + "/"

# Files to be uploaded
print("Copying files")
files = ['src/main.py', 'src/lib/core.py', 'src/lib/inout.py', 'src/lib/bounds.py']
for file in files:
    file_remote = remote_path + file
    file_local = local_path + file

    print(file_local + ' >>> ' + file_remote)
    try:
        sftp.remove(file_remote)
    except IOError:
        print(file + " was not on the cluster")
        pass

    sftp.put(file_local, file_remote)

print("Excecuting...")

# Create a local file that will be sent to the server (the infamous '.job' file)
# Main output folder:

time.sleep(.250)
current_time = time.time()
timestamp = datetime.datetime.fromtimestamp(current_time).strftime('%Y%m%d_%H:%M:%S')
current_folder = "run__" + timestamp

# È necessario creare la cartella /out/ in remoto!
try:
    sftp.mkdir(remote_path + "out/")
except IOError:
    # Se è già stata creata, non occorre ri-crearla
    pass

# Dato che la cartella corrente è un timestamp, siamo sicuri di poterla creare sempre (in remoto)
sftp.mkdir(remote_path + "out/" + current_folder)

with open("commands.job", "w") as fp:
    fp.write("#!/bin/bash \n")
    fp.write("#$ -q gpu \n")
    # fp.write("pip3 install --user networkx \n")
    # fp.write("pip3 install --user pandas \n")
    fp.write("source /nfsd/opt/anaconda3/anaconda.sh \n")
    fp.write("conda activate /nfsd/opt/anaconda3/tensorflow \n")


    # Formatting/constructing the instruction to be given:
    instruction = "time python3 " + remote_path + "src/main.py"
    # Options to be added:
    for k in parameters:
        if k != "prob" and k != "bound":
            instruction += " --" + k + " " + str(parameters[k])
    if parameters["prob"]:
        instruction += " --prob "
        if parameters["bound"]:
            instruction += " --bound"
    # current outputfolder for now is not used in our code
    # instruction += " --outfolder " + output_folder

    # Timeout is not currently used in our code
    # instruction += " --timeout 600 " + str(time_out)

    # Saving the output to a log file:
    output_logfilename = 'k=' + str(parameters['k']) + '_' + 'method=' + parameters['method']
    instruction += ' > '
    instruction += remote_path + "out/" + current_folder + '/' + output_logfilename + '_results.log'
    instruction += '\n'
    fp.write(instruction)

# Put the file on the current folder on the cluster and delete the local one
print(local_path + 'src/commands.job' ' >>> ' + remote_path + 'out/' + current_folder + '/commands.job')
sftp.put(local_path + 'src/commands.job', remote_path + 'out/' + current_folder + '/commands.job')
os.remove("commands.job")

# Give this job to the cluster
ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("export SGE_ROOT=/usr/share/gridengine \n" +  # necessary
                                                     "cd {0}out/{1} \n".format(remote_path,
                                                                               current_folder) +  # also necessary
                                                     "qsub -cwd commands.job")
# dev'essere "tutto assieme"  o si dimentica dell'export.
# una singola chiamata di exec_command fa dimenticare tutto quello che è stato fatto prima
# qsub -cwd == "current working directory". DEVE ESSERE MESSO PRIMA!!!

# Print output and errors
print(ssh_stdout.read().decode('utf-8'))
print(ssh_stderr.read().decode('utf-8'))

time.sleep(5 - .250)

sftp.close()
ssh.close()

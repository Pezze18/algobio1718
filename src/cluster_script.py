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
parameters['proteins_input']="../../../"+parameters['proteins_input']
parameters['samples_input']="../../../"+parameters['samples_input']
parameters['genes_input']="../../../"+parameters['genes_input']
parameters['filter_input']="../../../"+parameters['filter_input']
parameters['bestVectors']="../../../"+parameters['bestVectors']

server = "login.dei.unipd.it"

# Crea un'istanza del client SSH
ssh = paramiko.SSHClient()
# <<Required, since "login.dei.unipd.it" is not a "well known ssh host">> (???)
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Per evitare di scrivere il proprio username, se la macchina dalla quale lo sia avvia è "nota"
if getpass.getuser() == 'DavideDP':
    #pwd = str(getpass.getpass(prompt="Password: "))
    file=open("password.txt","r")
    pwd=file.readline()
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
files = ['src/main.py', 'src/lib/core.py', 'src/lib/inout.py', 'src/lib/bounds.py', 'src/lib/combinatorial.py', 'src/lib/auxiliary_functions.py']
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


def cluster_script(parameters):
    print("Excecuting...")
    ks = [7,8,9,10]

    waitingFor=open("../out/waitingFor.txt","a")
    waitingFor.write("\n")

    folder=remote_path + "out/"
    # È necessario creare la cartella /out/ in remoto!
    try:
        sftp.mkdir(folder)
    except IOError:        # Se è già stata creata, non occorre ri-crearla
        pass

    folder=folder+parameters["strategy"]
    try:
        sftp.mkdir(folder)
    except IOError:
        pass

    folder=folder + "/"+parameters["method"]
    try:
        sftp.mkdir(folder)
    except IOError:
        pass


    for k in ks:
        parameters['k']=k

        # Create a local file that will be sent to the server (the infamous '.job' file)
        time.sleep(.250)
        current_time = time.time()
        timestamp = datetime.datetime.fromtimestamp(current_time).strftime('%Y%m%d_%H:%M:%S')
        current_folder = folder+"/k="+str(parameters['k'])+"_run__" + timestamp

        # Dato che la cartella corrente è un timestamp, siamo sicuri di poterla creare sempre (in remoto)
        sftp.mkdir(current_folder )

        with open("commands.job", "w", newline='\n') as fp:
            fp.write("#!/bin/bash \n")
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

            # Saving the output to a log file:
            output_logfilename = 'k=' + str(parameters['k']) + '_' + 'method=' + parameters['method']
            instruction += ' > '
            instruction += current_folder + '/' + output_logfilename + '_results.log'
            instruction += '\n'
            fp.write(instruction)

        # Put the file on the current folder on the cluster and delete the local one
        print(local_path + 'src/commands.job' ' >>> ' + current_folder + '/commands.job')
        sftp.put(local_path + 'src/commands.job',  current_folder + '/commands.job')
        os.remove("commands.job")

        instruction=""
        with open("exec", "w", newline='\n') as exec:
            instruction+="export SGE_ROOT=/usr/share/gridengine \n"#non serve
            instruction+="cd {0} \n".format(current_folder)
            instruction+="qsub -q Q@runner-04 -cwd commands.job"
            exec.write(instruction)
            exec.close()
            with open("exec", "r", newline='\n') as exec:
                # Give this job to the cluster
                ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(exec.read())
        if os.path.exists("exec"):
            os.remove("exec")

        # Print output and errors
        line=ssh_stdout.read().decode('utf-8')
        start=line.find("b")
        end=line.find("(")
        idJob=line[start+1:end]
        k=parameters["k"]

        prob="det"
        if parameters["prob"]:
            prob="prob"

        line2=parameters["strategy"]+" "+prob+" "+parameters["method"]+" k="+str(k)+"_run__"+timestamp+" idJob:"+idJob+"\n"
        waitingFor.write(line2)
        print(line2)
        print(line)
        line=ssh_stderr.read().decode('utf-8')
        print(line)


    #close
    time.sleep(5 - .250)
    sftp.close()
    ssh.close()
    waitingFor.close()

if __name__ == "__main__":
    cluster_script(parameters)


# How to use
# 1. Login to https://app.globus.org/settings/developers and copy a project app id and secret
# 2. Use the id and secret to create and endpoint https://funcx.readthedocs.io/en/latest/sdk.html#client-credentials-with-clients
#     $ export FUNCX_SDK_CLIENT_ID="b0500dab-ebd4-430f-b962-0c85bd43bdbb"
#     $ export FUNCX_SDK_CLIENT_SECRET="ABCDEFGHIJKLMNOP0123456789="
# 3. Set up an endpoint on the computer that will run the tests, using these instructions: https://funcx.readthedocs.io/en/latest/endpoints.html
# 4. Create install-test.sh and run-test.sh on target computer

from globus_compute_sdk import Executor
import sys
import os

machine = sys.argv[1]
name = sys.argv[2]
branch = sys.argv[3]
endpoint = sys.argv[4]

with open(machine+'/env.sh', 'r') as file:
    env_file = file.read()

with open(machine+'/install.sh', 'r') as file:
    install_file = file.read()

with open(machine+'/run.sh', 'r') as file:
    run_file = file.read()

def run_on_endpoint(name, branch, env_file, install_file, run_file):
    import subprocess

    with open(name+"-test/env.sh", "w") as text_file:
        text_file.write("%s" % env_file)
        text_file.close()

    with open(name+"-test/install.sh", "w") as text_file:
        text_file.write("%s" % install_file)
        text_file.close()
    
    with open(name+"-test/run.sh", "w") as text_file:
        text_file.write("%s" % run_file)
        text_file.close()

    install_command = "cd {0}-test && chmod +x install.sh && ./install.sh {1} 2>&1 | tee build.log".format(name, branch)
    install_result = subprocess.run([install_command], shell=True, encoding="utf_8", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    if install_result.returncode != 0:
    return (install_result, None)

    run_command = "cd {0}-test && chmod +x run.sh && ./run.sh 2>&1 | tee test.log".format(name)
    run_result = subprocess.run([run_command], shell=True, encoding="utf_8", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    return (install_result, run_result)

gce = Executor(endpoint_id = endpoint)
future = gce.submit(run_on_endpoint, name, branch, env_file, install_file, run_file)
result = future.result()

with open("Build.log", "w") as text_file:
    text_file.write("%s" % result[0].stdout)
    text_file.close()

if result[1] != None:
    with open("Test.log", "w") as text_file:
        text_file.write("%s" % result[1].stdout)
        text_file.close()

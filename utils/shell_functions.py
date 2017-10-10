from subprocess import Popen, PIPE

def shell_true(string):
    p = Popen(string, stdout=PIPE, stderr=PIPE, shell=True)
    p.wait()
    stdout, stderr = p.communicate()
    print(stdout, stderr)

def shell_false(list):
    p = Popen(list, stdout=PIPE, stderr=PIPE)
    p.wait()
    stdout, stderr = p.communicate()
    print(stdout, stderr)
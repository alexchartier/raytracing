import subprocess
import platform
import datetime as dt
import pandas as pd

sTime = dt.datetime(2014, 5, 23, 6)
eTime = dt.datetime(2014, 5, 24)
freqMin = 2E6
freqMax = 30E6
distMin = 100
distMax = 10000
nLinks_req = 500
outFn = 'out.txt'

mac_executable = './wspr/WSPR.Rocks.Client.Tester/bin/Release/net6.0/osx.11.0-x64/publish/getWsprData'
linux_executable = './wspr/WSPR.Rocks.Client.Tester/bin/Release/net6.0/linux-x64/publish/getWsprData'

if 'mac' in platform.platform():
    executable = mac_executable
else:
    executable = linux_executable

cmd = '{0} "{1}" "{2}" "{3}" "{4}" {5} {6} {7} {8}'.format(
    executable, 
    sTime.strftime('%Y-%m-%d %H:%M:%S'), 
    eTime.strftime('%Y-%m-%d %H:%M:%S'),
    str(freqMin),
    str(freqMax), 
    str(distMin),
    str(distMax),
    nLinks_req,
    outFn,
)

print(cmd)
subprocess.call(cmd, shell=True)



wsprLinkDataAll = pd.read_csv('outputdftest.txt')


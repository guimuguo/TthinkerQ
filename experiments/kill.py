#!/usr/bin/env python
import subprocess
import os
ps = subprocess.Popen(['ipcs'], stdout=subprocess.PIPE).communicate()[0].split('\n')

mark = 0
idx = 0
while idx < len(ps):
    
    if  ps[idx].find('Message Queues') != -1:
        idx += 2
        mark = 1
        continue;

    if mark == 1:
        tokens = ps[idx].split()
        if len(tokens) == 6:
            id = tokens[1]
            os.system('ipcrm -q %s' % id)
            print(id, ' has been deleted.')

    idx += 1

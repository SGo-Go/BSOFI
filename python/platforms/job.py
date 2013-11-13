#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Module for running jobs on different platforms.
# ------------------------------------------------------------------

class Job(object):
    def __init__(self, cmd, **kw):
        super(Job, self).__init__()
        self.cmd = cmd
        self.attrs = kw
        pass

    def __getitem__(self, k):
        if k in self.attrs.keys():
            return self.attrs[k]
        return super(Job, self).__getitem__(k)

class JobManager(object):
    def __init__(self):
        super(JobManager, self).__init__()
        pass

    def submit(self, job):
        pass

class PBSJobManager(JobManager):
    def __init__(self, job_template):
        super(PBSJobManager, self).__init__()
        self.job_template = job_template
        pass

    def submit(self, job):
        from popen2 import popen2
        import time
        output, input = popen2('qsub  -S /bin/bash')
        job_template_handle = open(self.job_template, "r")
        job_script = job_template_handle.read()
        job_template_handle.close()
        job_script = job_script.replace(r"@cmd@", job.cmd)
        input.write(job_script)
        input.close()
        print "Job is submited: [%s]" % output.read()
        time.sleep(0.1)
        pass

class BashJobManager(JobManager):
    def __init__(self, cmd_submit = 'nohup'):
        super(BashJobManager, self).__init__()
        self.cmd_submit = cmd_submit
        pass

    def submit(self, job):
        from popen2 import popen2
        import time
        output, input = popen2(r'%s %s &'% (self.cmd_submit, job.cmd))
        #import os
        #os.system(r'%s %s &'% (self.cmd_submit, job.cmd))
        print "Job is submited: [%s]" % output.read()
        time.sleep(0.1)
        pass


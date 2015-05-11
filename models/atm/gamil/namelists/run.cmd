# @ shell = /bin/sh
#@ job_type=parallel
#@ notify_user = wab 
#@environment = COPY_ALL ; MP_TIMEOUT=1200 ; MALLOCMULTIHEAP=true; \
#                MP_SHARED_MEMORY=yes ; MBX_SIZE=160000000; MP_EUILIB=ip; MP_EAGER_LIMIT=0; \
#                XLSMPOPTS= "stack=86000000" ; datasize=unlimited ; MP_RMPOOL=1; \
#                MP_NODES = 4   ; MP_PROCS=30 ; MP_RESD=no ; MP_STDOUTMODE=ordered ; \
#                MP_EAGER_LIMIT=65536  ; MP_INFOLEVEL =6
#@ wall_clock_limit = 36000000
#@ class = small
#@ requirements =  (Arch == "R6000") && (OpSys == "AIX51")
#@ step_name = step1ipen0
#@ executable = /bin/poe
#@ arguments = ./camts < namelist
#@ output = run.$(jobid).$(stepid).out
#@ error =  run.$(jobid).$(stepid).err
#@ node_usage = not_shared
#@ total_tasks = 30 
#@ checkpoint = no
#@ resources = ConsumableCpus(1) 
#@ queue

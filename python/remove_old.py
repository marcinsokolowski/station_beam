import os
import sys
import time

def clean_old( directory="/tmp/station_beam" , max_age=432000 , do_remove=False ) : # 5 days 
   if not os.path.exists(directory) :
      print("INFO : directory %s does not exist -> nothing to clean up" % (directory))
      return False

   current_ux=time.time()
   i=0
   for x in os.walk(directory) :
   
#      print("%d\t" % i)
#      print(x)
      
      if i >= 1 :
         subdir=x[0]
         
         idx=subdir.find(directory)
         print("DEBUG : idx = %d" % (idx))
         
         if idx == 0 : # just safety check if directory path is not found -> path is different and should not be considered for removal
            creation_ux=os.path.getmtime(subdir)
            age_seconds = (current_ux - creation_ux)
            print("DEBUG : checking directory %s created at uxtime = %d vs. current_ux = %d -> age = %d [sec]" % (subdir,creation_ux,current_ux,age_seconds))
            if age_seconds >= max_age :
               print("\t ---> older than maximum age = %d seconds -> removing" % (max_age))

               cmd=("rm -fr %s" % subdir)
               print("\t ---> command = %s" % (cmd))
               if do_remove :
                  print("\tExecuting command: %s" % (cmd))
                  os.system(cmd)
                  
         else :
            print("ERROR : in code could not find directory name %s in string %s -> better not to remove something important !!!" % (directory,subdir))
         
         
      
      i=i+1

   return True         
   
if __name__ == "__main__":
   directory="/tmp/station_beam"
   if len(sys.argv) > 1:
      directory = sys.argv[1]
      
   max_age=432000
   if len(sys.argv) > 2:   
      max_age = float( sys.argv[2] )
      
   do_remove=False
   if len(sys.argv) > 3:      
      if int( sys.argv[3] ) > 0 :
         do_remove=True

   print("INFO : cleaning directory %s started at unixtime = %d , removing files and directories older than %d seconds (do_remove = %s)" % (directory,time.time(),max_age,do_remove))
   
   clean_old( directory=directory , max_age=max_age , do_remove=do_remove )
      
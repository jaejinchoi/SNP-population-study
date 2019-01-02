import os, sys, getopt


def show_help():

    print '-h, show_help'
    print '-S [save_foler_path]'
    #print '-l [float], ld begin'
    print '-L, save log, [SAVE_FOLDER]/log'
    print '-p [program_path]'
    print '-H [path pf hapmap]'
    print '-i, ignore kill code'
    print '-t [path], treemix index path'
    print '-R [int], max migration'
    print '-r [str], root-specific in treemix index'
    sys.exit()



#ld_value=ld_value_begin
try:
    opts, args  = getopt.getopt(sys.argv[1:], 'S:l:Lp:H:h:iLt:iR:r:P:')

except:
    show_help()


pre_save_folder_path=''
program_path=''
hapmap_path=''
prefix_path='./out'
treemix_index_path=''
treemix_migration_n=10
treemix_root='San'

kill_code=''
save_log_flag=False

for opt, arg in opts:

    if (opt=='-S'):
        pre_save_folder_path=os.path.abspath(arg)

    elif (opt=='-p'):
        program_path=os.path.abspath(arg)

    elif (opt=='-H'):
        hapmap_path=os.path.abspath(arg)

    elif (opt=='-h'):
        show_help()

    elif (opt=='-i'): #ignore kill code
        kill_code='-i'

    elif (opt=='-L'):
        save_log_flag=True

    elif (opt=='-t'):
        treemix_index_path=os.path.abspath(arg)

    elif (opt=='-R'):
        treemix_migration_n=int(arg)

    elif (opt=='-r'):
        treemix_root=str(arg)

    elif (opt=='-P'):
        prefix_path=str(arg)


ld_param=0
miss_param=0
maf_param=0
coverage_param=0 #read converagei, not necessary when using hapmap

'''
if (os.path.exists(save_folder_path)==False):
    os.makedirs(save_folder_path)
'''

#prefix_path=''
log_path=''

#3 parameters
for ld_param in [0.05, 0.5, 0]:#[0.2, 0.1, 0.05, 0.5, 0]: #[0, 0.5, 0.1, 0.2, 0.05]:

    for maf_param in [0.01, 0.05, 0.06]: #[0.05, 0.06, 0.1, 1]: #0.05 or 0.06 based on smallest clade or number of rare allele frequency of reports
        
        for miss_param in [0]: #[0.1, 0]: #wonder why miss_param make difference, even if there is no missing 'NN' case
            
            #prefix_path = '%g-%g-%g' % (ld_param, maf_param, miss_param)
            #prefix_path='out'
            save_folder_path = '%s/l%g-m%g' % (pre_save_folder_path, ld_param, maf_param) #need to separate for phylip(infile)
            fit_prefix_path='%s/%s' % (save_folder_path, prefix_path)

            if (save_log_flag==True):
                log_path= '> %s/%s' % (save_folder_path, 'log') #save-optional

            if (os.path.exists(save_folder_path)==False):
                os.makedirs(save_folder_path)

            #bootstrap on, -M and -p is 0
            print "cd %s && sh %s -l %g -m %g -p 0 -M 0 -P %s -b -H %s -t %s -R %d -r %s %s" % (save_folder_path, program_path, ld_param, maf_param, fit_prefix_path, hapmap_path, treemix_index_path, treemix_migration_n, treemix_root, log_path)
            

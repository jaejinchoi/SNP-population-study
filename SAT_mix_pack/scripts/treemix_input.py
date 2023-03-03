import sys, os, getopt


def clade_define(load_path='', individual_index_dict={}):
    
    taxon_clade_dict={}

    clade_name=''

    read_f = open(load_path, 'r')

    for read_line in read_f:
        read_line = read_line.strip()

        if (read_line!=''):
            if (read_line[0]=='/'):
                clade_name = read_line[1:]

                if (taxon_clade_dict.has_key(clade_name)==False):
                    taxon_clade_dict[clade_name]=[]

            else:
                taxon_clade_dict[clade_name].append(individual_index_dict[read_line])

    read_f.close()

    return taxon_clade_dict



def clade_define_mod(load_path=''):
    
    taxon_clade_dict={}

    str_list=[]

    clade_name=''

    read_f = open(load_path, 'r')

    for read_line in read_f:
        read_line = read_line.strip()

        if (read_line==''):
            continue

        else:
            str_list=read_line.split('\t')

            clade_name = str_list[1].replace(' ', '_')

            if (taxon_clade_dict.has_key(clade_name)==False):
                taxon_clade_dict[clade_name]=['%s_%s' % (str_list[2], str_list[4])]

            else:
                taxon_clade_dict[clade_name].append('%s_%s' % (str_list[2], str_list[4]))

    read_f.close()

    return taxon_clade_dict


def assign_individual_index(individual_list=[]):

    individual_index_dict={}

    index_cnt=0

    for n_item in individual_list:
        
        individual_index_dict[n_item]=index_cnt #should be 0 base
        index_cnt+=1

    return individual_index_dict



def count_population_allele(alleles_type_list=[], allele_distribute_list=[], population_class_dict={}):
    #alleles_type_list = 'ref/alt'

    pop_allele_str=''
    pop_allele_count_list=[]

    for n_population in population_class_dict:
        pop_allele_str=''
        pop_allele_count_list=[]

        for n_individual_index in population_class_dict[n_population]: #create single string of selected position

            pop_allele_str+=allele_distribute_list[n_individual_index] #0 base

        for n_allele_type in alleles_type_list:
            pop_allele_count_list.append(str(pop_allele_str.count(n_allele_type)))

        #fit for monoallele position
        if (len(pop_allele_count_list)==1):
            pop_allele_count_list.append('0')
        
        sys.stdout.write('%s ' % (','.join(pop_allele_count_list))) #delimiter is space

    sys.stdout.write('\n') #new line
        


def parse_hapmap(load_path='', index_path=''):

    read_str_list=[]

    #allele_consist_list=[]
    alleles_type_list=[] #biallelic

    individual_index_dict={}
    population_class_dict={}

    read_f = open(load_path, 'r')

    for read_line in read_f:
        read_line = read_line.strip()

        if (read_line==''):
            continue

        else:

            read_str_list = read_line.split('\t')

            if (read_line[:3]=='rs#'): #first row = header with individual id
                #individual_list = read_str_list[11:]
                individual_index_dict = assign_individual_index(read_str_list[11:])
                population_class_dict = clade_define(index_path, individual_index_dict)

                #print population_class_dict

                for n_pop in population_class_dict:
                    #print n_pop
                    sys.stdout.write('%s ' % (n_pop)) #delimiter is space

                sys.stdout.write('\n')
                
            else:
                alleles_type_list = read_str_list[1].split('/')

                count_population_allele(alleles_type_list, read_str_list[11:], population_class_dict)
                '''
                if (len(alleles_type_list)>1): #accept only biallelic case
                    #count_population_allele(alleles_type_list, read_str_list[11:], population_class_dict)
                    continue

                else:
                    print alleles_type_list
                '''


    read_f.close()

            

def show_help():
    print '-i [path], index_path'
    print '-h, show_help'
    sys.exit()

           
if __name__=='__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:hv')

    except:
        show_help()

    index_path=''

    for opt, arg in opts:

        if (opt=='-i'):
            index_path = os.path.abspath(arg)
            
        elif (opt=='-h'):
            show_help()

    ''' temporary use
    clade_dict = clade_define_mod(index_path)

    for n_pop in clade_dict:
        print '/%s' % (n_pop)

        for n_item in clade_dict[n_pop]:
            print n_item
    '''


    if (len(args)>0):

        for n_arg in args:
            load_path = os.path.abspath(n_arg)

            parse_hapmap(load_path, index_path)


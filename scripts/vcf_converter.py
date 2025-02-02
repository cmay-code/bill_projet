import os 
import re

parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
data_dir = str(parent_dir) + "/data/"



def convert_file(vcf_loc): #Converts the vcf input into a csv file that is stored in the data folder 
    sequence_dir = str(data_dir) + str(vcf_loc)

    file = open(sequence_dir, mode='r')
    vcfContents = file.readlines()[48:] #Ignore first 48 lines
    file.close()

    file_name = (re.findall("([^/]*$)", vcf_loc)[0]) 
    full_path_to_new_csv = (str(data_dir) + "files_converted_from_vcf_to_csv/" + str(file_name) + ".csv")

    with open(str(full_path_to_new_csv), "w") as file: #Overwrites the input file, in case the file already exists
        pass
    
    #Add " as a string delimeter in ordder to convert into csv
    for csvContent in vcfContents:
    
        csvContent = re.sub(r'^', '"', csvContent)

        csvContent = re.sub(r'\n', '', csvContent)
  
        csvContent = re.sub(r'\t$', '"', csvContent)

        if not re.search(r'"$', csvContent):
            csvContent = re.sub(r'$', '"', csvContent)

        #Turn all whitespaces to ","str(data_dir) + str(filename)
        csvContent = re.sub(r'\t', '","', csvContent)

        
        with open(full_path_to_new_csv, "a") as csvFile:
            csvFile.write(csvContent)
            csvFile.write('\n')
    
    csvFile.close()
    return ("csv file created in data folder")

print(convert_file("P15/KHV-U/P15-10.trimed1000.sv_sniffles.vcf"))


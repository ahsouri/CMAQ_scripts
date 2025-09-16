def read_csv(infiles,subsets,state):

    import numpy as np
    import csv
    from datetime import datetime
    
  
    output = {} ; keys = []
    first_file = True
    for infile in infiles:
        print(infile)
        with open(infile) as csv_file:
             csv_reader = csv.reader(csv_file, delimiter=',')
             counter = 0
             first = True
             for row in csv_reader:
                if(first):
                    for name in row:
                        if first_file:
                           output[name] = [] ; keys.append(name) ; nvar = len(keys)
                    first = False
                else:
                    is_in_state = False
                    for n in range(nvar):
                        if keys[n] == 'State Name':
                            if row[n] == state: is_in_state = True
                    if is_in_state:
                       for n in range(nvar):
                           if keys[n] in subsets:
                              if keys[n]=='Date Local':
                                 try:
                                     output[keys[n]].append(str(datetime.strptime(row[n], "%Y-%m-%d").strftime('%Y%m%d')))
                                 except:
                                     output[keys[n]].append(str(datetime.strptime(row[n], "%m/%d/%Y").strftime('%Y%m%d')))
                              elif keys[n]=='Time Local':
                                   output[keys[n]].append(str(datetime.strptime(row[n], "%H:%M").strftime("%H%M%S")))
                              else:       
                                   print(row[n])             
                                   output[keys[n]].append(row[n]) 
        first_file = False
    return output
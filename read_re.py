import re
import Levenshtein

with open('re_w_cutsites.txt','r') as REtable:
    re_dict = {}
    for line in REtable:
        line = line.rstrip()
        cleanline = line[:-1]
        cleanline_sliced = re.sub(r' (\S+)? +',r' ',cleanline)
        re_cutsite_pair = cleanline_sliced.split(' ')
        k = re_cutsite_pair[0]
        v = re_cutsite_pair[1]
        re_dict[k] = v
    # print(re_dict)

re_names = []
for key,val in re_dict.items():
    re_names.append(key)

input_re = input("Please enter your restriction enzyme of choice: ")
matches = [res_enz for res_enz in re_names if Levenshtein.distance(res_enz, input_re) < 2]

if not matches :
    print("Please enter a valid restriction enzyme. Check out a list of restriction enzymes here: https://rebase.neb.com/rebase/link_bionetc.")
elif input_re in matches:
    print(re_dict[input_re])
elif input_re not in matches:
    separator = ", "
    my_string = separator.join(matches)
    print(f'Do you mean {my_string}?')

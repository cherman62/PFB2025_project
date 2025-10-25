import re
import Levenshtein

def re_match(restriction_enzymes):
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
        re_names = []
        for key,val in re_dict.items():
            re_names.append(key)

    matches_list = []
    for res_enzyme in restriction_enzymes:
        matches = []
        for res_enz in re_names: 
            if Levenshtein.distance(res_enz, res_enzyme) < 2:
                matches.append(res_enz)
        matches_list.append(matches)

    results = []
    for count,res_enz_input in enumerate(restriction_enzymes):
        if not matches_list[count]:
            invalid_warning = 'Please enter a valid restriction enzyme. Check out a list of restriction enzymes here: https://rebase.neb.com/rebase/link_bionetc.'
            return invalid_warning
        elif res_enz_input in matches_list[count]:
            results.append(re_dict[res_enz_input])
        elif res_enz_input not in matches_list[count]:
            separator = ", "
            my_string = separator.join(matches_list[count])
            suggestion_warning = f'Do you mean {my_string}?'
            return suggestion_warning
    
    return results

restriction_enzymes = ['ZraI','BamHI','XmaI']
print(re_match(restriction_enzymes))

# def re_match(restriction_enzyme):
#     with open('re_w_cutsites.txt','r') as REtable:
#         re_dict = {}
#         for line in REtable:
#             line = line.rstrip()
#             cleanline = line[:-1]
#             cleanline_sliced = re.sub(r' (\S+)? +',r' ',cleanline)
#             re_cutsite_pair = cleanline_sliced.split(' ')
#             k = re_cutsite_pair[0]
#             v = re_cutsite_pair[1]
#             re_dict[k] = v
#         re_names = []
#         for key,val in re_dict.items():
#             re_names.append(key)
#     matches_list = []
#     for item in restriction_enzyme:
#         matches = res_enz for res_enz in re_names if Levenshtein.distance(res_enz, restriction_enzyme) < 3
#         print(matches)
#         matches_list.append(matches)
#     return matches_list
    
    # if not matches :
    #     invalid_warning = 'Please enter a valid restriction enzyme. Check out a list of restriction enzymes here: https://rebase.neb.com/rebase/link_bionetc.'
    #     return invalid_warning
    # elif restriction_enzyme in matches:
    #     return re_dict[restriction_enzyme]
    # elif restriction_enzyme not in matches:
    #     separator = ", "
    #     my_string = separator.join(matches)
    #     suggestion_warning = f'Do you mean {my_string}?'
    #     return suggestion_warning

# restrictionenz = ['BamHI','EcoRI']

# print(re_match(restrictionenz))
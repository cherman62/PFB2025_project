import re
import Levenshtein

def re_match(restriction_enzyme):
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
    matches = [res_enz for res_enz in re_names if Levenshtein.distance(res_enz, restriction_enzyme) < 3]
    if not matches :
        invalid_warning = 'Please enter a valid restriction enzyme. Check out a list of restriction enzymes here: https://rebase.neb.com/rebase/link_bionetc.'
        return invalid_warning
    elif restriction_enzyme in matches:
        return re_dict[restriction_enzyme]
    elif restriction_enzyme not in matches:
        separator = ", "
        my_string = separator.join(matches)
        suggestion_warning = f'Do you mean {my_string}?'
        return suggestion_warning

restrictionenz = 'SmaI'
print(re_match(restrictionenz))

# def re_dict_names():
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
#         return re_dict, re_names

# re_dictionary_namelist = re_dict_names()
# #This variable is a tuple. The first element is a dictionary {'Enzyme': 'Cutsite', ...} and the second element is a list of res enzymes ['BamHI','EcoRI', ...]
# print(re_dictionary_namelist[0])
# print('')
# print('')
# print(re_dictionary_namelist[1])
# print('')
# print('')

# def re_match(restriction_enzyme):
#     matches = [res_enz for res_enz in re_dictionary_namelist[1] if Levenshtein.distance(res_enz, restriction_enzyme) < 3]
#     if not matches :
#         invalid_warning = 'Please enter a valid restriction enzyme. Check out a list of restriction enzymes here: https://rebase.neb.com/rebase/link_bionetc.'
#         return invalid_warning
#     elif restriction_enzyme in matches:
#         return re_dictionary_namelist[0][restriction_enzyme]
#     elif restriction_enzyme not in matches:
#         separator = ", "
#         my_string = separator.join(matches)
#         suggestion_warning = f'Do you mean {my_string}?'
#         return suggestion_warning

# restrictionenz = 'BamHSJKHSFKI'
# print(re_match(restrictionenz))
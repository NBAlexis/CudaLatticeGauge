from MesonStructures import all_signs, all_deltas

all_signs_lst = all_signs()
all_delta_lst = all_deltas()

for i in range(20):
    all_type = len(all_signs_lst[i])
    for j in range(all_type):
        print(f"===={i}-{j}====")
        print(all_signs_lst[i][j])
        print(all_delta_lst[i][j])


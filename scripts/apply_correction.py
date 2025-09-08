from apply_lacs_correction import apply_selected_offsets_and_note
import json
import numpy as np

def read_json_result(fname):
    with open(fname, "r", encoding="utf-8") as f:
        report2 = json.load(f)
    print(len(report2))
    for bmrb_id in report2:
        if len(report2[bmrb_id])==0:
            print(bmrb_id)
        for list_id in report2[bmrb_id]:
            if len(report2[bmrb_id][list_id])==0:
                print(bmrb_id, list_id)
            ca = []
            cb = []
            c = []
            n = []
            ha = []
            ofs={}
            for method in report2[bmrb_id][list_id]:
                for atom in report2[bmrb_id][list_id][method]:
                    if report2[bmrb_id][list_id][method][atom] is not None:
                        if atom=='ca': ca.append(report2[bmrb_id][list_id][method][atom])
                        if atom=='cb': cb.append(report2[bmrb_id][list_id][method][atom])
                        if atom=='c': c.append(report2[bmrb_id][list_id][method][atom])
                        if atom=='n': n.append(report2[bmrb_id][list_id][method][atom])
                        if atom=='ha':ha.append(report2[bmrb_id][list_id][method][atom])

            if np.median(ca)>0.2:
                #print (bmrb_id,list_id,ca,np.mean(ca),np.std(ca),'CA')
                ofs['CA']=round(np.median(ca),2)
            if np.median(cb)>0.2:
                #print (bmrb_id,list_id,cb,np.mean(cb),np.std(cb),'CB')
                ofs['CB']=round(np.median(cb),2)
            if np.median(c)>0.2:
                #print (bmrb_id,list_id,c,np.mean(c),np.std(c),'C')
                ofs['C']=round(np.median(c),2)
            if np.median(n)>0.2:
                #print (bmrb_id,list_id,n,np.mean(n),np.std(n),'N')
                ofs['N']=round(np.median(n),2)
            if ofs:
                if 'CA' in ofs and 'CB' not in ofs:
                    avg= round((ofs['CA']+np.median(cb))/2.0,2)
                    ofs['CA']=avg
                    ofs['CB']=avg
                if 'CB' in ofs and 'CA' not in ofs:
                    avg = round((ofs['CB'] + np.median(ca)) / 2.0, 2)
                    ofs['CA'] = avg
                    ofs['CB'] = avg
                parts = [f"{a}={ofs[a]:+g}" for a in ofs]
                details = (
                        f"CS reference correction applied to list_id {list_id} : "
                        + (", ".join(parts) if parts else "none")
                        + ". Source: PyLACS."
                )
                print (details)
                flg = [True if i>5.0 else False for i in ofs.values()]
                if True in flg:
                    out_path=f'../scratch/corrected/bmr{bmrb_id}_corrected_wa.str'
                else:
                    out_path=f'../scratch/corrected/bmr{bmrb_id}_corrected.str'
                counts = apply_selected_offsets_and_note(
                    input_path=f'/reboxitory/2025/08/BMRB/macromolecules/bmr{bmrb_id}/bmr{bmrb_id}_3.str',
                    output_path=out_path,
                    list_id=int(list_id),
                    offsets=ofs,
                    atoms=list(ofs.keys()),

                    release_author='BMRB',
                    release_details=details,
                )
                print (counts)
if __name__ == '__main__':
    read_json_result('bmrb_offsets.json')

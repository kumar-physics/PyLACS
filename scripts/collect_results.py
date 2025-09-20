import os
import json
import pandas as pd
import plotly.express as px


def read_json_in_folders(root_dir):
    """
    Reads JSON files inside every folder under the given root directory.

    Parameters:
        root_dir (str): Path to the root directory
    """
    before={}
    after={}
    id=[]
    atms=['ca','cb','c','n']
    for atom in atms:
        before[atom]=[]
        after[atom]=[]
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if file.lower().endswith(".json"):
                file_path = os.path.join(dirpath, file)
                try:
                    with open(file_path, "r") as f:
                        (head, tail) = os.path.split(file_path)
                        corrected_file = f'{os.path.split(head)[0]}/corrected/{tail.split("_")[0]}/{tail}'
                        if os.path.exists(corrected_file):
                            data = json.load(f)
                            print(f"\n--- Contents of {os.path.basename(file_path)} ---")
                            (head,tail)=os.path.split(file_path)
                            print (head,tail,os.path.split(head))
                            #json_data = json.dumps(data)
                            #print(json.dumps(data, indent=4))
                            for list_id in data:
                                id.append(tail.split("_")[0])
                                for atom in atms:
                                    try:
                                        before[atom].append(data[list_id]['offsets'][atom])
                                    except KeyError:
                                        before[atom].append(None)
                            corrected_file = f'{os.path.split(head)[0]}/corrected/{tail.split("_")[0]}/{tail}'
                            print (corrected_file)
                            with open(corrected_file,'r') as f2:
                                data2 = json.load(f2)
                                for list_id in data2:
                                    for atom in atms:
                                        try:
                                            after[atom].append(data2[list_id]['offsets'][atom])
                                        except KeyError:
                                            after[atom].append(None)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")

    print (len(id))
    for atom in atms:
        print (len(before[atom]),len(after[atom]))
    df = pd.DataFrame({
        "BMRB_ID" : id,
        "CA(before)":before['ca'],
        "CA(after)":after['ca'],
        "CB(before)": before['cb'],
        "CB(after)": after['cb'],
        "C(before)": before['c'],
        "C(after)": after['c'],
        "N(before)": before['n'],
        "N(after)": after['n']
    })
    fig = px.scatter(df,x='CA(before)',y="CA(after)",hover_name="BMRB_ID",marginal_x="histogram",marginal_y="histogram")
    fig.show()
    fig.write_html('ca_correlation.html')
    fig = px.scatter(df, x='CB(before)', y="CB(after)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('cb_correlation.html')
    fig = px.scatter(df, x='C(before)', y="C(after)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('c_correlation.html')
    fig = px.scatter(df, x='N(before)', y="N(after)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('n_correlation.html')
    fig = px.scatter(df, x='CA(before)', y="CB(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('ca_cb_correlation.html')
    fig = px.scatter(df, x='CA(before)', y="C(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('ca_c_correlation.html')
    fig = px.scatter(df, x='CA(before)', y="N(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('ca_n_correlation.html')
    fig = px.scatter(df, x='C(before)', y="N(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('c_n_correlation.html')

# Example usage:
if __name__ == "__main__":
    root_directory = "/home/nmrbox/kbaskaran/lacs/bayes"
    read_json_in_folders(root_directory)

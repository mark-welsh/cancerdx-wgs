import argparse
import json
from ordered_set import OrderedSet

parser = argparse.ArgumentParser()
parser.add_argument("--top_panel", dest="top_panel", required=True)
parser.add_argument("--mid_panel", dest="mid_panel", required=True)
parser.add_argument("--main_id", dest="main_id", required=True, help='"main" ID (e.g. tumor sample, DGD ID)')
parser.add_argument("--secondary_id", dest="secondary_id", required=True, help='"secondary" ID (e.g. normal sample, case ID)')
args = parser.parse_args()

with open(args.top_panel, 'r') as top_in:
    top_lines = top_in.readlines()[1:]

with open(args.mid_panel, 'r') as bot_in:
    bot_lines = bot_in.readlines()[1:]

plot = {
    "cnv_data": [],
    "cnv_layout": [],
    "baf_data": [],
    "baf_layout": [],
}

chr_labels = OrderedSet()

x = list()
y = list()
color = list()
text = list() 

baf_x = list()
baf_y = list()
baf_color = list()
baf_text = list()

chrm_bp = list()
start_chrm = None
for i, line in enumerate(top_lines):
    cols = line.strip().split('\t')
    chrm = cols[0]
    chr_labels.add('chr'+ chrm if 'chr' not in chrm else chrm)
    if start_chrm is None:
        start_chrm = chrm
    if chrm != start_chrm:
        chrm_bp.append(i)
        start_chrm = chrm
    log2 = cols[3]
    if cols[4] == 'normal':
        c = 'gray'
    elif cols[4] == 'gain':
        c = 'darkgreen'
    elif cols[4] == 'loss':
        c = 'red'
    else:
        raise TypeError
    gene = cols[5]
    last_chrm = chrm

    x.append(i)
    y.append(log2)
    color.append(c)
    text.append(gene)

plot["cnv_data"].append({
    "x": x,
    "y": y,
    "marker": {
        "color": color,
        "size": "5"
    },
    "text": text,
    "type": "scattergl",
    "mode": "markers",
    "hoverinfo": "text"
})

plot["cnv_layout"] = {
    "legend": {"orientation": 'h', "bordercolor": '#D3D3D3', "borderwidth": 1},
    "margin": {"l": 50, "r": 50, "b": 50, "t": 50},
    "yaxis": {"range": [-3.0, 3.0], "zeroline": False, "title": 'log2ratio'},
    "xaxis": {"range": [-100, plot["cnv_data"][0]["x"][-1]], "showline": False, "showgrid": False, "tickmode": "array",
        "ticktext": list(chr_labels), "showticklabels": True, "tickfont": {"size": 10}
    },
    "hovermode": 'closest',
    "shapes": [
        {
            "type": "line",
            "y0": 0.4, "y1": 0.4,
            "x0": 0, "x1": plot["cnv_data"][0]["x"][-1],
            "line": {"color": "darkgreen"}
        },
        {
            "type": "line",
            "y0": -0.7, "y1": -0.7,
            "x0": 0, "x1": plot["cnv_data"][0]["x"][-1],
            "line": {"color": "red"}
        }
    ]
}

for i, line in enumerate(bot_lines):
    cols = line.strip().split('\t')
    baf_x.append(cols[0])
    baf_y.append(cols[1])
    baf_text.append(cols[3])
    if cols[4] == 'normal':
        c = 'gray'
    elif cols[4] == 'gain':
        c = 'darkgreen'
    elif cols[4] == 'loss':
        c = 'red'
    else:
        raise TypeError
    baf_color.append(c)

plot["baf_data"].append({
    "x": baf_x,
    "y": baf_y,
    "marker": {
        "color": baf_color,
        "size": "5"
    },
    "text": baf_text,
    "type": "scattergl",
    "mode": "markers",
    "hoverinfo": "text"
})

plot["baf_layout"] = {
    "yaxis": {"range": [-0.025, 1.025], "zeroline": False, "title": 'B-allele frequency'},
    "xaxis": {
        "range": [-100, plot["cnv_data"][0]["x"][-1]], "showline": False, "showgrid": False, "tickmode": "array",
        "ticktext": list(chr_labels), "showticklabels": True, "tickfont": {"size": 10}
    },
    "hovermode": 'closest',
    "shapes": []
}

chrm_ticks = list() 
for i, bp in enumerate(chrm_bp):
    if i == 0:
        chrm_ticks.append(int(bp/2))
    elif i+1 == len(chrm_bp):
        chrm_ticks.append(int((chrm_bp[-1]+bp)/2))
    else:
        chrm_ticks.append( int((bp - chrm_bp[i-1])/2) + chrm_bp[i-1])
    chrm_line = {
        "type": "line",
        "y0": -3.0, "y1": 3.0,
        "x0": int(bp), "x1": int(bp),
        "line": {"color": "#D3D3D3"}
    }
    plot["cnv_layout"]["shapes"].append(chrm_line)
    plot["baf_layout"]["shapes"].append(chrm_line)


plot["cnv_layout"]["xaxis"]["tickvals"] = chrm_ticks
plot["baf_layout"]["xaxis"]["tickvals"] = chrm_ticks

plot["dgd_id"] = args.main_id
plot["case_id"] = args.secondary_id

with open('{}_{}.json'.format(args.main_id, args.secondary_id), 'w') as fout:
    json.dump(plot, fout, indent=2)

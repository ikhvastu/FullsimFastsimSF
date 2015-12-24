from ROOT import kRed, kCyan, kGreen
from ROOT import TGraph
from ROOT import TMultiGraph
#from ROOT import TGraphErrors
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import TLatex
from array import array
import numpy as np

folder = "results_ptmore25/"
file = open(folder + 'felm.txt', 'r')
good_value = []
x = []
y = []
ex = []
ey = []
prev_value = float(0.)
for i in file:
    value = i.split(" ")
    if prev_value != float(value[1]):
        x.append(float(value[1]))
        ex.append(float(value[2]))
        y.append(float(value[3]))
        ey.append(float(value[4]))
        prev_value = float(value[1])
    else:
        continue

n = len(x)
x_a = array("d", x)
y_a = array("d", y)
ex_a = array("d", ex)
ey_a = array("d", ey)
gr = TGraph(n, y_a, x_a)
file.close()

file = open(folder + 'feli.txt', 'r')
good_value = []
x = []
y = []
ex = []
ey = []
prev_value = float(0.)
for i in file:
    value = i.split(" ")
    if prev_value != float(value[1]):
        x.append(float(value[1]))
        ex.append(float(value[2]))
        y.append(float(value[3]))
        ey.append(float(value[4]))
        prev_value = float(value[1])
    else:
        continue

n = len(x)
x_a1 = array("d", x)
y_a1 = array("d", y)
ex_a1 = array("d", ex)
ey_a1 = array("d", ey)
gr1 = TGraph(n, y_a1, x_a1)

file = open(folder + 'fmum.txt', 'r')
good_value = []
x = []
y = []
ex = []
ey = []
prev_value = float(0.)
for i in file:
    value = i.split(" ")
    if prev_value != float(value[1]):
        x.append(float(value[1]))
        ex.append(float(value[2]))
        y.append(float(value[3]))
        ey.append(float(value[4]))
        prev_value = float(value[1])
    else:
        continue

n = len(x)
x_a = array("d", x)
y_a = array("d", y)
ex_a = array("d", ex)
ey_a = array("d", ey)
gr2 = TGraph(n, y_a, x_a)
file.close()

file = open(folder + 'fmui.txt', 'r')

x = []
y = []
ex = []
ey = []
prev_value = float(0.)
for i in file:
    value = i.split(" ")
    if prev_value != float(value[1]):
        x.append(float(value[1]))
        ex.append(float(value[2]))
        y.append(float(value[3]))
        ey.append(float(value[4]))
        prev_value = float(value[1])
    else:
        continue

n = len(x)
x_a1 = array("d", x)
y_a1 = array("d", y)
ex_a1 = array("d", ex)
ey_a1 = array("d", ey)

gr3 = TGraph(n, y_a1, x_a1)

file = open(folder + 'felmulti.txt', 'r')
x = []
y = []
ex = []
ey = []
for i in file:
    value = i.split(" ")
    x.append(float(value[1]))
    ex.append(float(value[2]))
    y.append(float(value[3]))
    ey.append(float(value[4]))

n = len(x)
x_a10 = np.array(float(x[0]))
y_a10 = np.array(float(y[0]))
gr0multi = TGraph(1, y_a10, x_a10)
gr0multi.SetMarkerColor(kCyan - 7)
gr0multi.SetMarkerStyle(21)

x_a11 = np.array(float(x[1]))
y_a11 = np.array(float(y[1]))
gr1multi = TGraph(1, y_a11, x_a11)
gr1multi.SetMarkerColor(kGreen - 4)
gr1multi.SetMarkerStyle(21)

x_a12 = np.array(float(x[2]))
y_a12 = np.array(float(y[2]))
gr2multi = TGraph(1, y_a12, x_a12)
gr2multi.SetMarkerColor(kRed - 6)
gr2multi.SetMarkerStyle(21)

file.close()


file = open(folder + 'fmumulti.txt', 'r')
x = []
y = []
ex = []
ey = []
for i in file:
    value = i.split(" ")
    x.append(float(value[1]))
    ex.append(float(value[2]))
    y.append(float(value[3]))
    ey.append(float(value[4]))

n = len(x)
x_a20 = np.array(float(x[0]))
y_a20 = np.array(float(y[0]))
gr3multi = TGraph(1, y_a20, x_a20)
gr3multi.SetMarkerColor(kCyan - 7)
gr3multi.SetMarkerStyle(21)

x_a21 = np.array(float(x[1]))
y_a21 = np.array(float(y[1]))
gr4multi = TGraph(1, y_a21, x_a21)
gr4multi.SetMarkerColor(kGreen - 4)
gr4multi.SetMarkerStyle(21)

x_a22 = np.array(float(x[2]))
y_a22 = np.array(float(y[2]))
gr5multi = TGraph(1, y_a22, x_a22)
gr5multi.SetMarkerColor(kRed - 6)
gr5multi.SetMarkerStyle(21)

file.close()

c1 = TCanvas("c1", "c1")
c2 = TCanvas("c2", "c2")
c3 = TCanvas("c3", "c3")
c1.Divide(2, 1)
c2 = c1.cd(1)
c2.SetLogx()

mg = TMultiGraph()

gr1.SetTitle("TGraphErrors Example")
gr1.SetLineColor(2)
gr1.SetMarkerStyle(1)

gr.SetMarkerColor(1)
gr.SetMarkerStyle(1)

mg.Add(gr1)
mg.Add(gr)
mg.Add(gr0multi)
mg.Add(gr1multi)
mg.Add(gr2multi)
mg.Draw("APL")
mg.GetXaxis().SetTitle("Background eff")
mg.GetXaxis().SetLimits(0.005, 1.0)
mg.GetXaxis().SetMoreLogLabels()
mg.GetXaxis().SetTitleSize(0.04)
mg.GetXaxis().SetLabelSize(0.04)
mg.GetXaxis().SetNoExponent()
mg.GetYaxis().SetTitle("Eff Signal(T1t^{4} 1.2/0.8)")
mg.GetYaxis().SetTitleOffset(1.3)
mg.GetYaxis().SetTitleSize(0.04)
mg.GetYaxis().SetLabelSize(0.04)
mg.SetMaximum(1.)
mg.SetMinimum(.5)

leg = TLegend(0.15, 0.75, 0.55, 0.9)
leg.SetTextSize(0.04)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(gr, "miniisolation", "l")
leg.AddEntry(gr1, "standard isolation", "l")
leg.AddEntry(gr0multi, "WP-L", "p")
leg.AddEntry(gr1multi, "WP-M", "p")
leg.AddEntry(gr2multi, "WP-T", "p")
leg.Draw("same")

latex = TLatex()
latex.SetTextSize(0.04)
#latex.DrawLatex(0.073, 0.93, "el, 10 GeV < p_T < 25 GeV")
latex.DrawLatex(0.073, 0.93, "el, p_T > 25 GeV")

c3 = c1.cd(2)
c3.SetLogx()

mg1 = TMultiGraph()

gr3.SetTitle("TGraphErrors Example")
gr3.SetLineColor(2)
gr3.SetMarkerStyle(1)

gr2.SetMarkerColor(1)
gr2.SetMarkerStyle(1)

mg1.Add(gr2)
mg1.Add(gr3)
mg1.Add(gr3multi)
mg1.Add(gr4multi)
mg1.Add(gr5multi)
mg1.Draw("APL")
mg1.GetXaxis().SetTitle("Background eff")
mg1.GetXaxis().SetLimits(0.005, 1.0)
mg1.GetXaxis().SetTitleSize(0.04)
mg1.GetXaxis().SetLabelSize(0.04)
mg1.GetXaxis().SetNoExponent()
mg1.GetXaxis().SetMoreLogLabels()
mg1.GetYaxis().SetTitle("Eff Signal(T1t^{4} 1.2/0.8)")
mg1.GetYaxis().SetTitleOffset(1.3)
mg1.GetYaxis().SetTitleSize(0.04)
mg1.GetYaxis().SetLabelSize(0.04)
mg1.SetMaximum(1.)
mg1.SetMinimum(.5)

leg1 = TLegend(0.15, 0.75, 0.55, 0.9)
leg1.SetTextSize(0.04)
leg1.SetFillStyle(0)
leg1.SetBorderSize(0)
leg1.AddEntry(gr2, "miniisolation", "l")
leg1.AddEntry(gr3, "standard isolation", "l")
leg1.AddEntry(gr3multi, "WP-L", "p")
leg1.AddEntry(gr4multi, "WP-M", "p")
leg1.AddEntry(gr5multi, "WP-T", "p")
leg1.Draw("same")

latex1 = TLatex()
latex1.SetTextSize(0.04)
#latex1.DrawLatex(0.055, 0.93, "mu, 10 GeV < p_T < 25 GeV")
latex1.DrawLatex(0.055, 0.93, "mu, p_T > 25 GeV")

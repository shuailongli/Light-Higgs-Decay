#include "ToolFunctions.h"

double abs_Delta_K(double x)
{
double limit = 0;
if ((0<=x)&&(x<=0.0393259)) {
      limit =0.233363+(x-0)*(0.233321-0.233363)/(0.0393259-0);
}
if ((0.0393259<=x)&&(x<=0.0763876)) {
      limit =0.233321+(x-0.0393259)*(0.235169-0.233321)/(0.0763876-0.0393259);
}
if ((0.0763876<=x)&&(x<=0.113448)) {
      limit =0.235169+(x-0.0763876)*(0.235127-0.235169)/(0.113448-0.0763876);
}
if ((0.113448<=x)&&(x<=0.143097)) {
      limit =0.235127+(x-0.113448)*(0.235093-0.235127)/(0.143097-0.113448);
}
if ((0.143097<=x)&&(x<=0.18016)) {
      limit =0.235093+(x-0.143097)*(0.240812-0.235093)/(0.18016-0.143097);
}
if ((0.18016<=x)&&(x<=0.255344)) {
      limit =0.240812+(x-0.18016)*(0.246624-0.240812)/(0.255344-0.18016);
}
if ((0.255344<=x)&&(x<=0.320998)) {
      limit =0.246624+(x-0.255344)*(0.254634-0.246624)/(0.320998-0.255344);
}
if ((0.320998<=x)&&(x<=0.387716)) {
      limit =0.254634+(x-0.320998)*(0.27153-0.254634)/(0.387716-0.320998);
}
if ((0.387716<=x)&&(x<=0.459727)) {
      limit =0.27153+(x-0.387716)*(0.287211-0.27153)/(0.459727-0.387716);
}
if ((0.459727<=x)&&(x<=0.530681)) {
      limit =0.287211+(x-0.459727)*(0.306262-0.287211)/(0.530681-0.459727);
}
if ((0.530681<=x)&&(x<=0.600577)) {
      limit =0.306262+(x-0.530681)*(0.331892-0.306262)/(0.600577-0.530681);
}
if ((0.600577<=x)&&(x<=0.6673)) {
      limit =0.331892+(x-0.600577)*(0.365525-0.331892)/(0.6673-0.600577);
}
if ((0.6673<=x)&&(x<=0.70649)) {
      limit =0.365525+(x-0.6673)*(0.399384-0.365525)/(0.70649-0.6673);
}
if ((0.70649<=x)&&(x<=0.738265)) {
      limit =0.399384+(x-0.70649)*(0.425956-0.399384)/(0.738265-0.70649);
}
if ((0.738265<=x)&&(x<=0.774276)) {
      limit =0.425956+(x-0.738265)*(0.454287-0.425956)/(0.774276-0.738265);
}
if ((0.774276<=x)&&(x<=0.808174)) {
      limit =0.454287+(x-0.774276)*(0.504458-0.454287)/(0.808174-0.774276);
}
if ((0.808174<=x)&&(x<=0.839953)) {
      limit =0.504458+(x-0.808174)*(0.551206-0.504458)/(0.839953-0.808174);
}
if ((0.839953<=x)&&(x<=0.867502)) {
      limit =0.551206+(x-0.839953)*(0.632182-0.551206)/(0.867502-0.839953);
}
if ((0.867502<=x)&&(x<=0.90035)) {
      limit =0.632182+(x-0.867502)*(0.742805-0.632182)/(0.90035-0.867502);
}
if ((0.90035<=x)&&(x<=0.926843)) {
      limit =0.742805+(x-0.90035)*(0.865796-0.742805)/(0.926843-0.90035);
}
if ((0.926843<=x)&&(x<=0.948048)) {
      limit =0.865796+(x-0.926843)*(1.05925-0.865796)/(0.948048-0.926843);
}
if ((0.948048<=x)&&(x<=0.974559)) {
      limit =1.05925+(x-0.948048)*(1.40482-1.05925)/(0.974559-0.948048);
}
if ((0.974559<=x)&&(x<=0.992532)) {
      limit =1.40482+(x-0.974559)*(1.14804-1.40482)/(0.992532-0.974559);
}
if ((0.992532<=x)&&(x<=0.995673)) {
      limit =1.14804+(x-0.992532)*(0.886719-1.14804)/(0.995673-0.992532);
}
if ((0.995673<=x)&&(x<=1.01364)) {
      limit =0.886719+(x-0.995673)*(0.695979-0.886719)/(1.01364-0.995673);
}
if ((1.01364<=x)&&(x<=1.03797)) {
      limit =0.695979+(x-1.01364)*(0.573354-0.695979)/(1.03797-1.01364);
}
if ((1.03797<=x)&&(x<=1.06653)) {
      limit =0.573354+(x-1.03797)*(0.4839-0.573354)/(1.06653-1.03797);
}
if ((1.06653<=x)&&(x<=1.09934)) {
      limit =0.4839+(x-1.06653)*(0.435633-0.4839)/(1.09934-1.06653);
}
if ((1.09934<=x)&&(x<=1.1311)) {
      limit =0.435633+(x-1.09934)*(0.39536-0.435633)/(1.1311-1.09934);
}
if ((1.1311<=x)&&(x<=1.16285)) {
      limit =0.39536+(x-1.1311)*(0.370583-0.39536)/(1.16285-1.1311);
}
if ((1.16285<=x)&&(x<=1.19674)) {
      limit =0.370583+(x-1.16285)*(0.367544-0.370583)/(1.19674-1.16285);
}
if ((1.19674<=x)&&(x<=1.23379)) {
      limit =0.367544+(x-1.19674)*(0.347293-0.367544)/(1.23379-1.19674);
}
if ((1.23379<=x)&&(x<=1.26873)) {
      limit =0.347293+(x-1.23379)*(0.350048-0.347293)/(1.26873-1.23379);
}
if ((1.26873<=x)&&(x<=1.30897)) {
      limit =0.350048+(x-1.26873)*(0.358557-0.350048)/(1.30897-1.26873);
}
if ((1.30897<=x)&&(x<=1.3365)) {
      limit =0.358557+(x-1.30897)*(0.349933-0.358557)/(1.3365-1.30897);
}
if ((1.3365<=x)&&(x<=1.37039)) {
      limit =0.349933+(x-1.3365)*(0.352711-0.349933)/(1.37039-1.3365);
}
if ((1.37039<=x)&&(x<=1.41062)) {
      limit =0.352711+(x-1.37039)*(0.344206-0.352711)/(1.41062-1.37039);
}
if ((1.41062<=x)&&(x<=1.44133)) {
      limit =0.344206+(x-1.41062)*(0.338645-0.344206)/(1.44133-1.41062);
}
if ((1.44133<=x)&&(x<=1.47944)) {
      limit =0.338645+(x-1.44133)*(0.330483-0.338645)/(1.47944-1.44133);
}
if ((1.47944<=x)&&(x<=1.51226)) {
      limit =0.330483+(x-1.47944)*(0.317361-0.330483)/(1.51226-1.47944);
}
if ((1.51226<=x)&&(x<=1.54931)) {
      limit =0.317361+(x-1.51226)*(0.299875-0.317361)/(1.54931-1.51226);
}
if ((1.54931<=x)&&(x<=1.58531)) {
      limit =0.299875+(x-1.54931)*(0.28565-0.299875)/(1.58531-1.54931);
}
if ((1.58531<=x)&&(x<=1.61284)) {
      limit =0.28565+(x-1.58531)*(0.276539-0.28565)/(1.61284-1.58531);
}
if ((1.61284<=x)&&(x<=1.65307)) {
      limit =0.276539+(x-1.61284)*(0.263415-0.276539)/(1.65307-1.61284);
}
if ((1.65307<=x)&&(x<=1.68906)) {
      limit =0.263415+(x-1.65307)*(0.255003-0.263415)/(1.68906-1.65307);
}
if ((1.68906<=x)&&(x<=1.71871)) {
      limit =0.255003+(x-1.68906)*(0.246867-0.255003)/(1.71871-1.68906);
}
if ((1.71871<=x)&&(x<=1.75682)) {
      limit =0.246867+(x-1.71871)*(0.235154-0.246867)/(1.75682-1.71871);
}
if ((1.75682<=x)&&(x<=1.79176)) {
      limit =0.235154+(x-1.75682)*(0.227645-0.235154)/(1.79176-1.75682);
}
if ((1.79176<=x)&&(x<=1.82881)) {
      limit =0.227645+(x-1.79176)*(0.220374-0.227645)/(1.82881-1.79176);
}
if ((1.82881<=x)&&(x<=1.89869)) {
      limit =0.220374+(x-1.82881)*(0.198357-0.220374)/(1.89869-1.82881);
}
if ((1.89869<=x)&&(x<=1.96115)) {
      limit =0.198357+(x-1.89869)*(0.187404-0.198357)/(1.96115-1.89869);
}
if ((1.96115<=x)&&(x<=2)) {
      limit =0.187404+(x-1.96115)*(0.177079-0.187404)/(2-1.96115);
}
return limit;
}
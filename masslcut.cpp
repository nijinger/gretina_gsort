TCutG *load_masslcut(){
//========= Macro generated from object: masslcut/Graph
//========= by ROOT version5.34/36
   
   TCutG *cutg = new TCutG("masslcut",18);
   cutg->SetVarX("fthetaL[0]");
   cutg->SetVarY("dT[0]");
   cutg->SetTitle("Graph");
   cutg->SetFillColor(1);
   cutg->SetPoint(0,0.323236,-420.43);
   cutg->SetPoint(1,0.304877,-306.989);
   cutg->SetPoint(2,0.295697,-197.911);
   cutg->SetPoint(3,0.379461,-106.285);
   cutg->SetPoint(4,0.595181,-35.7475);
   cutg->SetPoint(5,0.925645,-3.75118);
   cutg->SetPoint(6,1.21824,39.8801);
   cutg->SetPoint(7,1.33987,61.6958);
   cutg->SetPoint(8,1.42708,76.2396);
   cutg->SetPoint(9,1.44085,51.5152);
   cutg->SetPoint(10,1.34676,32.6083);
   cutg->SetPoint(11,1.13907,-5.93275);
   cutg->SetPoint(12,0.950889,-43.7466);
   cutg->SetPoint(13,0.786804,-155.006);
   cutg->SetPoint(14,0.681239,-272.084);
   cutg->SetPoint(15,0.600918,-390.616);
   cutg->SetPoint(16,0.433391,-421.885);
   cutg->SetPoint(17,0.323236,-420.43);
   return cutg;
}

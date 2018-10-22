TCutG *load_massrcut(){
//========= Macro generated from object: massrcut/Graph
//========= by ROOT version5.34/36
   
   auto cutg = new TCutG("massrcut",15);
   cutg->SetVarX("fthetaR[0]");
   cutg->SetVarY("dT[0]");
   cutg->SetTitle("Graph");
   cutg->SetFillColor(1);
   cutg->SetPoint(0,0.298623,436.656);
   cutg->SetPoint(1,0.279339,264.914);
   cutg->SetPoint(2,0.293113,174.273);
   cutg->SetPoint(3,0.502479,64.5494);
   cutg->SetPoint(4,0.835813,12.0729);
   cutg->SetPoint(5,1.19118,-30.8624);
   cutg->SetPoint(6,1.42259,-61.0762);
   cutg->SetPoint(7,1.42534,-30.8624);
   cutg->SetPoint(8,1.13609,16.8435);
   cutg->SetPoint(9,0.868871,99.5338);
   cutg->SetPoint(10,0.744904,215.618);
   cutg->SetPoint(11,0.648485,342.834);
   cutg->SetPoint(12,0.574105,427.114);
   cutg->SetPoint(13,0.406061,466.869);
   cutg->SetPoint(14,0.298623,436.656);
   return cutg;
}

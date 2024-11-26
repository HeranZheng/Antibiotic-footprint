clear all
clc

cd('C:\Users\zheng\OneDrive - University College London\Research\3#Application\2024\Antibiotic');

%% readdata
Full = readtable('PCU_intensity_1', 'FileType', 'spreadsheet', 'Sheet', 'result');
Bridge = readtable('Bridge', 'FileType', 'spreadsheet', 'Sheet', 'Sheet1');

Continetal = unique(Full.Continetal);
Class = unique(Full.Class);
Type=unique(Full.Item);

Full_selected_R={};
Full_selected_R1={};

for k=1:size(Class,1)
       Full_selected = Full(strcmp(Full.Class,Class(k)),:);
       Bridge_selected = Bridge(strcmp(Bridge.Type,Class(k)),:);

    for i=1:size(Continetal,1)

        
%% Proxy-benchmark
        Proxy=Bridge_selected(strcmp(Bridge_selected.Country,Continetal(i)),3);
            if string(Proxy.Proxy) == 'Full' %% Directly use

                        Full_selected_f=Full_selected(strcmp(Full_selected.Continetal,Continetal(i)),:);

                        Full_selected_R=cat(1,Full_selected_R,Full_selected_f);

                                                         
            else
                    if string(Proxy.Proxy) == 'Null' 
                            

                             else
                           Full_selected_2=Full(strcmp(Full.Area,Proxy.Proxy),:);% finding the benchmark for the continental
                           Share_1=table2array(Full_selected_2(:,26:end))./(Full_selected_2.Anti_forecast); % column 26 must have a number

                           Benckmark.(Type{1})=Share_1(1:11,:) ;               
                           Benckmark.(Type{2})=Share_1(11+1:11+11,:) ;   
                           Benckmark.(Type{3})=Share_1(11+11+1:11+11+11,:) ;  
%%        
                            Full_selected_1=Full_selected(strcmp(Full_selected.Continetal,Continetal(i)),:);

                            Full_selected_12=Full_selected_1(~isnan(Full_selected_1.Tetracyclines_kg),:);% sample

                            Full_selected_11=Full_selected_1(isnan(Full_selected_1.Tetracyclines_kg),:);

                            List=unique(Full_selected_11.Area);
                                    
                                 if size(List,1)==0
                                      %% No countries for the category
                                  else

                                     for j=1:size(List,1)
                                        Share=[];
                                         Temp_s={};
                                
                                          Full_selected_f = Full_selected_11(strcmp(Full_selected_11.Area,List(j,1)),:);
                                           Type_Item = unique(Full_selected_f.Item);
                                   
                                              for tt1=1:size(Type_Item,1)
                                                   Temp_s(tt1,1) =  Type(strcmp(Type,Type_Item(tt1)),:);
                                               end

                                                  for t=1:size(Temp_s,1)

                                                    Share((t-1)*11+1:(t-1)*11+11,:)= Benckmark.(string(Temp_s(t)));

                                                    end 
                                                        Full_selected_f{:,26:end}=Full_selected_f{:,24}.*Share; % crucial
        
                                                      Full_selected_R=cat(1,Full_selected_R,Full_selected_f);
                                          
                                        end

                                

                                end

                                        Full_selected_R1=cat(1,Full_selected_R1,Full_selected_12);% merge estimated data and sample data

                            
            end
    end
    
    end

end

Full_selected_R2=cat(1,Full_selected_R,Full_selected_R1);

%% Sort out known value

%Full_selected_Intermediate=Full_selected_R2(find(cell2mat(table2cell(Full_selected_R2(:,25)))~=0),:);

%Full_selected=cell2mat(table2cell(Full_selected_Intermediate(:,24))).*(cell2mat(table2cell(Full_selected_Intermediate(:,26:end)))./cell2mat(table2cell(Full_selected_Intermediate(:,25))));

%Full_selected_R2{find(cell2mat(table2cell(Full_selected_R2(:,25)))~=0),26:41}=(Full_selected);

%% preparation
Regions_Fabio = readtable('regions', 'FileType', 'auto');
Regions_EXBASE_BRIDGE = readtable('regions', 'FileType', 'auto','Sheet','Bridge');
region=readtable('Bridge','FileType', 'auto','Sheet','Country_EX');
Regions_EXBASE_BRIDGE=cell2mat(table2cell(Regions_EXBASE_BRIDGE(:,2:end)));
Regions_EXBASE_BRIDGE(isnan(Regions_EXBASE_BRIDGE))=0;

Year=[2010:2020];
Inventory=struct()

for i = 1:size(Year,2)

    Anti=Full_selected_R2(find(Full_selected_R2.YearCode==Year(i)),:);
    
    for j=1:size(Type,1)

         Anti_1 = Anti(strcmp(Anti.Item,Type(j)),:);
    
         for k=1 :size(Regions_Fabio,1)
                if isempty(Anti_1(strcmp(Anti_1.FabioGroup,Regions_Fabio.iso3c(k)),24))

            Inventory_1(k,1:17)=repmat(mat2cell(zeros(1),1,1),1,17);

                else

            Inventory_1(k,1)=Anti_1(strcmp(Anti_1.FabioGroup,Regions_Fabio.iso3c(k)),24);

            Inventory_1(k,2:17)=Anti_1(strcmp(Anti_1.FabioGroup,Regions_Fabio.iso3c(k)),26:41);

            Name_in(k,1)=Regions_Fabio.iso3c(k);  % not in use                

                 end
         end
        Inventory_2=cell2mat(table2cell(Inventory_1));
        Inventory_2(isnan(Inventory_2))=0;
        Inventory_FB.(['Y_',num2str(Year(i))]).([Type{j}])=Inventory_2;% Convert to FABIO
        Inventory_ex.(['Y_',num2str(Year(i))]).([Type{j}])=Regions_EXBASE_BRIDGE'*Inventory_2;% Convert to EXIOBASE
    end
end

save('Anti_Inventory','Inventory_ex');

%% Linking with EXIOBASE
cd('C:\Users\zheng\OneDrive - University College London\Research\3#Application\2024\Antibiotic');
%% deflation
Deflator=readtable('./Data/Deflator/CPI.csv','FileType', 'auto');
for k=1:49
  % Country Benchmark     
     Constant_F1=Deflator(strcmp(Deflator.DataSource,region{k,1}),:);
      
     deflator((k-1)*200+1:(k-1)*200+200,:)=repmat(table2array(Constant_F1(:,5:15))./100,200,1);
end

Type_country=readtable('Bridge','FileType', 'auto','Sheet','Type2');
%Type_country=readtable('Bridge','FileType', 'auto','Sheet','Type'); Gloria
Sector_country=readtable('Bridge','FileType', 'auto','Sheet','Sector_bridge');
Bridge_sector=cell2mat(table2cell(Sector_country(:,2:end)));
Bridge_sector(isnan(Bridge_sector))=0;

cd('C:\Users\zheng\OneDrive - University College London\Research\#Data\IO Data\MRIO\EXIOBASE');

old=pwd;


for i=10:size(Year,2)

path=['.\IOT_' char(string(Year(i))) '_pxp'];
cd(path)

X=importdata("x.txt");
X1=X.data./deflator(:,i);
X2=X.textdata;

Y=importdata("Y.txt");
Y1=Y.data./deflator(:,i);
Y2=Y.textdata;

Z=importdata("Z.txt");
Z1=Z.data./deflator(:,i);
Z2=Z.textdata;

A1=Z1./X1';
A1(isnan(A1))=0;

%A=importdata("A.txt");
%A1=A.data;
%A2=A.textdata;

I=eye(size(A1));

L=inv(I-A1);

%% Antiobiotic inventory
n_country=49;
n_sector=200;
E=zeros(17,n_sector*n_country);
% cattle-96, pig-100, poutry-101
IN_1=(Inventory_ex.(['Y_' char(string(Year(1,i)))]).Cattle)';
IN_2=(Inventory_ex.(['Y_' char(string(Year(1,i)))]).Pig)';
IN_3=(Inventory_ex.(['Y_' char(string(Year(1,i)))]).Chicken)';

IN_1(isnan(IN_1))=0;
IN_2(isnan(IN_2))=0;
IN_3(isnan(IN_3))=0;

for j=1:n_country
E(:,(j-1)*n_sector+9)=IN_1(:,j);% Cattle
E(:,(j-1)*n_sector+10)=IN_2(:,j);% Pig
E(:,(j-1)*n_sector+11)=IN_3(:,j);% Chicken
end

E_I=E./X1';
E_I(isnan(E_I))=0;

E_I1=E_I(1,:);

for y=1:49
E_I2(y,:)=E_I1(:,(y-1)*200+1:(y-1)*200+200);
end

E3((i-1)*49+1:(i-1)*49+49,:)=mean(E_I2,2);



multi=diag(E_I(1,:))*L;
multi_1=(E_I(1,:))*L;

for i1=1:n_country
Footprint(:,i1)=multi*sum(Y1(:,(i1-1)*7+1:(i1-1)*7+7),2);

Footprint_1(i1,:)=multi_1*diag(sum(Y1(:,(i1-1)*7+1:(i1-1)*7+7),2));

% Full
Footprint_F=multi*diag(sum(Y1(:,(i1-1)*7+1:(i1-1)*7+7),2));

Footprint_F1=zeros(9800,200);
    for i2=1:n_country
        Footprint_F1=Footprint_F1+Footprint_F(:,(i2-1)*200+1:(i2-1)*200+200);
    end

Footprint_F2(:,(i1-1)*200+1:(i1-1)*200+200)=Footprint_F1;

end
%i1=31
%F_check(:,i)=Y1((i1-1)*200+1:(i1-1)*200+200,(i1-1)*7+1);

Outcome.(['F_',num2str(Year(i))]).('Full')=Footprint_F2;
Outcome.(['F_',num2str(Year(i))]).('Forward')=Footprint;
Outcome.(['F_',num2str(Year(i))]).('Backward')=Footprint_1;
Outcome.(['F_',num2str(Year(i))]).('Intenity')=mean(E_I2,2);
SDA.(['F_',num2str(Year(i))]).('L')=L;
SDA.(['F_',num2str(Year(i))]).('E')=E_I;
SDA.(['F_',num2str(Year(i))]).('Y')=Y1;

clearvars Footprint Footprint_1 Footprint_F2
cd(old)
end


cd('C:\Users\zheng\OneDrive - University College London\Research\3#Application\2024\Antibiotic');

%%% Inventory by type
for i=1:size(Year,2)

%% Antiobiotic inventory- exiobase
n_country=49;
n_sector=200;
E=zeros(17,n_sector*n_country);
% cattle-96, pig-100, poutry-101
IN_1=(Inventory_ex.(['Y_' char(string(Year(1,i)))]).Cattle)';
IN_2=(Inventory_ex.(['Y_' char(string(Year(1,i)))]).Pig)';
IN_3=(Inventory_ex.(['Y_' char(string(Year(1,i)))]).Chicken)';

IN_1(isnan(IN_1))=0;
IN_2(isnan(IN_2))=0;
IN_3(isnan(IN_3))=0;

for j=1:n_country
E(:,(j-1)*n_sector+9)=IN_1(:,j);% Cattle
E(:,(j-1)*n_sector+10)=IN_2(:,j);% Pig
E(:,(j-1)*n_sector+11)=IN_3(:,j);% Chicken
end

E_Share=E(2:end,:)./sum(E(2:end,:),1);
E_Share(isnan(E_Share))=0;
E_Share(isinf(E_Share))=0;
E_1{i}=E_Share;
%writematrix(E,['./Intensity/' char(string(Year(1,i))) '.csv']);

end


%%Result template
Result(:,1)=repmat(region,size(Year,2),1);
Result(:,2)=repmat(Type_country,size(Year,2),1);
for i=1:size(Year,2)
Result{(i-1)*49+1:(i-1)*49+49,3}=repmat(Year(i),49,1);
end


%% aggregate-forward linkage
for i=1:size(Year,2)
Footprint=Outcome.(['F_',num2str(Year(i))]).('Forward');

% aggregate
for i1=1:n_country
Footprint1(i1,:) = sum(Footprint((i1-1)*n_sector+1:(i1-1)*n_sector+n_sector,:),1);
end

Outcome.(['F_',num2str(Year(i))]).('Matrix49x49')=Footprint1;

for y=1:49
Footprint2(y,1)=Footprint1(y,y);% local
end

Result{(i-1)*49+1:(i-1)*49+49,4}=sum(Footprint1,2);% production
Result{(i-1)*49+1:(i-1)*49+49,5}=sum(Footprint1,1)';% consumption

Result{(i-1)*49+1:(i-1)*49+49,6}=Footprint2;% local
Result{(i-1)*49+1:(i-1)*49+49,7}=sum(Footprint1,1)'-Footprint2;% import-related
end

for x=1:49
    Outcome.F_2010.Matrix49x49(x,x)=0;
Outcome.F_2020.Matrix49x49(x,x)=0;
end

(Outcome.F_2010.Matrix49x49-(Outcome.F_2010.Matrix49x49)')./1000;
(Outcome.F_2020.Matrix49x49-(Outcome.F_2020.Matrix49x49)')./1000;





%%Result product template
Sector_10=readtable('Bridge','FileType', 'auto','Sheet','Sector');
Anti_type=readtable('Bridge','FileType', 'auto','Sheet','Anti');
Anti_Bridge=readtable('Bridge','FileType', 'auto','Sheet','Anti_Bridge');
Anti_type_short=Anti_Bridge.Properties.VariableNames;
Anti_type_short=Anti_type_short(:,2:7);

Anti_Bridge=cell2mat(table2cell(Anti_Bridge(:,2:end)));
Anti_Bridge(isnan(Anti_Bridge))=0;

n_product=size(Bridge_sector,2);
n_anti=size(Anti_type,1);

Result_product=table();
% country/region
for i=1:49
Result_product{(i-1)*n_product.*n_anti+1:(i-1)*n_product.*n_anti+n_product.*n_anti,1}=repmat(region{i,1},size(Bridge_sector,2).*size(Anti_type,1),1);
end
% region_type
for i=1:49
Result_product{(i-1)*n_product.*n_anti+1:(i-1)*n_product.*n_anti+n_product.*n_anti,2}=repmat(Type_country{i,1},size(Bridge_sector,2).*size(Anti_type,1),1);
end
Result_product=repmat(Result_product(1:n_product*49*n_anti,:),11,1);
% sector
for i=1:n_product
Result_product_t((i-1)*n_anti+1:(i-1)*n_anti+n_anti,1)=repmat(Sector_10(i,1),n_anti,1);
end
Result_product(:,3)=repmat(repmat(Result_product_t,49,1),size(Year,2),1);
% Year
for i=1:size(Year,2)
Result_product{(i-1)*n_product*49*n_anti+1:(i-1)*n_product*49*n_anti+n_product*49*n_anti,4}=repmat(Year(i),n_product*49*n_anti,1);
end
% Anti
Result_product{:,5}=repmat(table2cell(Anti_type),n_product*49*size(Year,2),1);


%% Layout for Sankey trade
Result_sankey=table();
Result_sankey_1=table();
% country/region Exporter
for i=1:49
Result_sankey{(i-1)*n_product+1:(i-1)*n_product+n_product,1}=repmat(region{i,1},size(Bridge_sector,2),1);
end
Result_sankey=repmat(Result_sankey,490,1);

% region_type_Exporter
for i=1:49
Result_sankey_1{(i-1)*n_product+1:(i-1)*n_product+n_product,1}=repmat(Type_country{i,1},size(Bridge_sector,2),1);
end
Result_sankey{:,2}=repmat(table2cell(Result_sankey_1),490,1);

% country/region Importer
for i=1:49
Result_sankey{(i-1)*490*10+1:(i-1)*490*10+490*10,3}=repmat(append(region{i,1},'1'),490*10,1);
end

% region_type_Importer
for i=1:49
Result_sankey{(i-1)*490*10+1:(i-1)*490*10+490*10,4}=repmat(append(Type_country{i,1},'1'),490*10,1);
end

% sector export
Result_sankey{:,5}=repmat(repmat(cellstr(table2cell(Sector_10)),49,1),490,1);

% sector import
sector_1=repmat(repmat(append(cellstr(table2cell(Sector_10))','1'),1,49),490,1);

Result_sankey{:,6}=sector_1(:);

% Anti
for i=1:size(Anti_type_short,2)
Result_sankey_Anti(:,i)=repmat(cellstr(Anti_type_short{1,i}),size(Result_sankey,1),1);
end
% Link
Result_sankey_total(:,1:6)=repmat(Result_sankey,size(Anti_type_short,2),1);

Result_sankey_total(:,7)=Result_sankey_Anti(:);


%% Pop
Pop=xlsread('Population','Population' );
Pop(isnan(Pop))=0;
Pop=Pop(2:end,:);
Bridge_pop=xlsread('regions','Bridge' );
Bridge_pop(isnan(Bridge_pop))=0;
Pop_1=Bridge_pop'*Pop;


%% aggregate-backward linkage 
for i= 1:size(Year,2)
%% Pop
for j=1:49
Pop_2((j-1)*n_product.*n_anti+1:(j-1)*n_product.*n_anti+n_product.*n_anti,i)=repmat(Pop_1(j,i),n_product.*n_anti,1);
end
    
%%Anti Total
Footprint_0=Outcome.(['F_',num2str(Year(i))]).('Full');

E_Share=E_1{i}';


for j=1:n_anti
Footprint_1=E_Share(:,j).*Footprint_0;

Footprint2=zeros(200,9800);

for i1=1:n_country
Footprint2=Footprint2+Footprint_1((i1-1)*n_sector+1:(i1-1)*n_sector+n_sector,:);
end

Footprint2=sum(Footprint2,1);


for i1=1:n_country
Footprint4(j,(i1-1)*n_product+1:(i1-1)*n_product+n_product)=(Footprint2(:,(i1-1)*n_sector+1:(i1-1)*n_sector+n_sector)*Bridge_sector);
end

sum(Footprint4,1);


%% Local
Footprint2=zeros(200,9800);

for i1=1:n_country
Footprint2(1:200,(i1-1)*200+1:(i1-1)*200+200) =Footprint_1((i1-1)*n_sector+1:(i1-1)*n_sector+n_sector,(i1-1)*n_sector+1:(i1-1)*n_sector+n_sector);
end

Footprint2=sum(Footprint2,1);

for i1=1:n_country
Footprint5(j,(i1-1)*n_product+1:(i1-1)*n_product+n_product)=(Footprint2(:,(i1-1)*n_sector+1:(i1-1)*n_sector+n_sector)*Bridge_sector);
end

%% Trade-imported
Footprint_2=Footprint_1;
Footprint2=zeros(200,9800);

for i2=1:n_country
    Footprint_2((i2-1)*200+1:(i2-1)*200+200,(i2-1)*200+1:(i2-1)*200+200)=0;
end

for i1=1:n_country
Footprint2 = Footprint2+Footprint_2((i1-1)*n_sector+1:(i1-1)*n_sector+n_sector,:);
end

Footprint2=sum(Footprint2,1);

for i1=1:n_country
Footprint6(j,(i1-1)*n_product+1:(i1-1)*n_product+n_product)=(Footprint2(:,(i1-1)*n_sector+1:(i1-1)*n_sector+n_sector)*Bridge_sector);
end


%% Trade-exported
Footprint_2=Footprint_1;
Footprint2=zeros(9800,200);
for i2=1:n_country
    Footprint_2((i2-1)*200+1:(i2-1)*200+200,(i2-1)*200+1:(i2-1)*200+200)=0;
end

for i1=1:n_country
Footprint2 = Footprint2+Footprint_2(:,(i1-1)*n_sector+1:(i1-1)*n_sector+n_sector);
end

Footprint2=sum(Footprint2,2);

for i1=1:n_country
Footprint7(j,(i1-1)*n_product+1:(i1-1)*n_product+n_product)=(Bridge_sector'*Footprint2((i1-1)*n_sector+1:(i1-1)*n_sector+n_sector,:))';
end

%%% trade flow for sankey
Footprint_2=Footprint_1;
for i2=1:n_country
    Footprint_2((i2-1)*200+1:(i2-1)*200+200,(i2-1)*200+1:(i2-1)*200+200)=0;
end

for i1=1:n_country
    for j1=1:n_country
    Footprint_sankey = Footprint_2(:,(i1-1)*n_sector+1:(i1-1)*n_sector+n_sector)*Bridge_sector;

    Footprint_sankey_1((j1-1)*10+1:(j1-1)*10+10,(i1-1)*10+1:(i1-1)*10+10)=Bridge_sector'* Footprint_sankey((j1-1)*n_sector+1:(j1-1)*n_sector+n_sector,:);
    end
end

Footprint_sankey_2(:,j)=Footprint_sankey_1(:);

end

Result_product{(i-1)*n_product*49*n_anti+1:(i-1)*n_product*49*n_anti+n_product*49*n_anti,6}=Footprint4(:);% total footprint
Result_product{(i-1)*n_product*49*n_anti+1:(i-1)*n_product*49*n_anti+n_product*49*n_anti,7}=Footprint5(:);% local
Result_product{(i-1)*n_product*49*n_anti+1:(i-1)*n_product*49*n_anti+n_product*49*n_anti,8}=Footprint6(:);% import
Result_product{(i-1)*n_product*49*n_anti+1:(i-1)*n_product*49*n_anti+n_product*49*n_anti,9}=Footprint7(:);% exported

end
%% Sankey
Footprint_sankey_2=Footprint_sankey_2*Anti_Bridge;

Result_sankey_total{:,8}=Footprint_sankey_2(:);

Result_sankey_total(find(table2array(Result_sankey_total(:,8))==0),:)=[];


% change name
Result_sankey_total.Properties.VariableNames(1) = "Exporter";
Result_sankey_total.Properties.VariableNames(2) = "Exporter_region";
Result_sankey_total.Properties.VariableNames(3) = "Importer";
Result_sankey_total.Properties.VariableNames(4) = "Importer_region";
Result_sankey_total.Properties.VariableNames(5) = "Export_product";
Result_sankey_total.Properties.VariableNames(6) = "Import_product";
Result_sankey_total.Properties.VariableNames(7) = "Type";
Result_sankey_total.Properties.VariableNames(8) = "Value";

Result_product{:,10}=Pop_2(:);

writetable(Result_product,'Result_product.xlsx','Sheet',1);
writetable(Result_sankey_total,'Result_product_Anti.csv');
writetable(Result,'Result.xlsx','Sheet',1);


%% SDA
Year_SDA=[2010,2020];
POP_SDA(:,1)=Pop_1(:,1);
POP_SDA(:,2)=Pop_1(:,7);
POP_SDA(:,3)=Pop_1(:,11);

for i=1:size(Year_SDA,2)-1

Delta_E=Outcome.(['F_',num2str(Year_SDA(i+1))]).('Forward')-Outcome.(['F_',num2str(Year_SDA(i))]).('Forward'); 
% E  

E_D=SDA.(['F_',num2str(Year_SDA(i+1))]).('E');
E_0=SDA.(['F_',num2str(Year_SDA(i))]).('E');

E_t=E_D(1,:)';
E_0=E_0(1,:)';

Y11=SDA.(['F_',num2str(Year_SDA(i+1))]).('Y');
Y00=SDA.(['F_',num2str(Year_SDA(i))]).('Y');

for j=1:49
    Y3(:,j)=sum(Y11(:,(j-1)*7+1:(j-1)*7+7),2);
    Y0(:,j)=sum(Y00(:,(j-1)*7+1:(j-1)*7+7),2);
end

D_E=E_t-E_0;

D_E1=diag(D_E)*SDA.(['F_',num2str(Year_SDA(i+1))]).('L')*Y3;
D_E2=diag(D_E)*SDA.(['F_',num2str(Year_SDA(i))]).('L')*Y0;

E_D=(D_E1+D_E2)/2;

for j=1:49
    E_D_1(j,:)=sum(E_D((j-1)*n_sector+1:(j-1)*n_sector+n_sector,:),1);
end

E_SDA=sum(sum(E_D,1));

%%%% L
L_E=SDA.(['F_',num2str(Year_SDA(i+1))]).('L')-SDA.(['F_',num2str(Year_SDA(i))]).('L'); 
    
      L_E1=diag(E_0)*L_E*Y3;
      L_E2=diag(E_t)*L_E*Y0;

L_D=(L_E1+L_E2)/2;

for j=1:49
    L_D_1(j,:)=sum(L_D((j-1)*n_sector+1:(j-1)*n_sector+n_sector,:),1);
end
L_SDA=sum(sum(L_D,1));

%%% F_S
FV_15=sum(sum(Y3,2));
FV_12=sum(sum(Y0,2));
%%Delta F_S

FS_15=Y3./FV_15;
FS_12=Y0./FV_12;

FS_D=FS_15-FS_12;
    
      FS_E1=diag(E_0)*SDA.(['F_',num2str(Year_SDA(i))]).('L')*(FS_D.*FV_15);
      FS_E2=diag(E_t)*SDA.(['F_',num2str(Year_SDA(i+1))]).('L')*(FS_D.*FV_12);

FS_D=(FS_E1+FS_E2)/2;

% 200 FS_D

FS_D11=zeros(200,49);

for j=1:49
    FS_D11=FS_D11+FS_D((j-1)*200+1:(j-1)*200+200,:);
end


for j=1:49
    FS_D_1(j,:)=sum(FS_D((j-1)*n_sector+1:(j-1)*n_sector+n_sector,:),1);
end

FS_SDA=sum(sum(FS_D,1));

%%% FV Per 
pop_2015=POP_SDA(:,i+1);
pop_2012=POP_SDA(:,i);

FVP_15=FV_15./pop_2015;
FVP_12=FV_12./pop_2012;

FVP_D=FVP_15-FVP_12;
    
      FV_E1=diag(E_0)*SDA.(['F_',num2str(Year_SDA(i))]).('L')*(FS_12.*(FVP_D.*pop_2015)');
      FV_E2=diag(E_t)*SDA.(['F_',num2str(Year_SDA(i+1))]).('L')*(FS_15.*(FVP_D.*pop_2012)');

FVP_D=(FV_E1+FV_E2)/2;

for j=1:49
    FVP_D_1(j,:)=sum(FVP_D((j-1)*n_sector+1:(j-1)*n_sector+n_sector,:),1);
end

FVP_SDA=sum(sum(FVP_D,1));

%%% Pop
pp_D=pop_2015-pop_2012;
    
      pp_E1=diag(E_0)*SDA.(['F_',num2str(Year_SDA(i))]).('L')*(FS_12.*(FVP_12.*pp_D)');
      pp_E2=diag(E_t)*SDA.(['F_',num2str(Year_SDA(i+1))]).('L')*(FS_15.*(FVP_15.*pp_D)');

pp_D=(pp_E1+pp_E2)/2;

for j=1:49
    pp_D_1(j,:)=sum(pp_D((j-1)*n_sector+1:(j-1)*n_sector+n_sector,:),1);
end

pp_SDA=sum(sum(pp_D,1));

SDA_detail.(['F_',num2str(Year_SDA(i))]).('Intensity')=E_D_1;
SDA_detail.(['F_',num2str(Year_SDA(i))]).('Trade')=L_D_1;
SDA_detail.(['F_',num2str(Year_SDA(i))]).('Finalstructure')=FS_D_1;
SDA_detail.(['F_',num2str(Year_SDA(i))]).('demand')=FVP_D_1;
SDA_detail.(['F_',num2str(Year_SDA(i))]).('pop')=pp_D_1;

SDA_Outcome(1,i)=E_SDA; % Intensity
SDA_Outcome(2,i)=L_SDA; % Leontif
SDA_Outcome(3,i)=FS_SDA; %structure
SDA_Outcome(4,i)=FVP_SDA;% per capita
SDA_Outcome(5,i)=pp_SDA;% pop

end

%% Driving factor
Driving_type=readtable('Bridge','FileType', 'auto','Sheet','factor');
n_driving=size(Driving_type,1);

Result_driving=table();
% region
for i=1:49
Result_driving{(i-1)*n_driving+1:(i-1)*n_driving+n_driving,1}=repmat(region{i,1},n_driving,1);
end
% region_type
for i=1:49
Result_driving{(i-1)*n_driving+1:(i-1)*n_driving+n_driving,2}=repmat(Type_country{i,1},n_driving,1);
end
% driving type
for i=1:49
Result_driving{(i-1)*n_driving+1:(i-1)*n_driving+n_driving,3}=table2cell(Driving_type);
end
% cat
Result_driving1=Result_driving;
for i=1:size(Year_SDA,2)-2
    
Result_driving=[Result_driving;Result_driving1];
end

% period
for i=1:size(Year_SDA,2)-1
Result_driving{(i-1)*n_driving*49+1:(i-1)*n_driving*49+n_driving*49,4}=repmat(i,49*5,1);
end


for i=1:size(Year_SDA,2)-1
Temp_E=SDA_detail.(['F_',num2str(Year_SDA(i))]).('Intensity');
Temp_L=SDA_detail.(['F_',num2str(Year_SDA(i))]).('Trade');
Temp_S=SDA_detail.(['F_',num2str(Year_SDA(i))]).('Finalstructure');
Temp_D=SDA_detail.(['F_',num2str(Year_SDA(i))]).('demand');
Temp_P=SDA_detail.(['F_',num2str(Year_SDA(i))]).('pop');


%% total
SDA_detail_total(1,:)=sum(SDA_detail.(['F_',num2str(Year_SDA(i))]).('Intensity'),1);
SDA_detail_total(2,:)=sum(SDA_detail.(['F_',num2str(Year_SDA(i))]).('Trade'),1);
SDA_detail_total(3,:)=sum(SDA_detail.(['F_',num2str(Year_SDA(i))]).('Finalstructure'),1);
SDA_detail_total(4,:)=sum(SDA_detail.(['F_',num2str(Year_SDA(i))]).('demand'),1);
SDA_detail_total(5,:)=sum(SDA_detail.(['F_',num2str(Year_SDA(i))]).('pop'),1);

Result_driving{(i-1)*n_driving*49+1:(i-1)*n_driving*49+n_driving*49,5}=SDA_detail_total(:);

%% Import
for j=1:49
    Temp_E(j,j)=0;
    Temp_L(j,j)=0;
    Temp_S(j,j)=0;
    Temp_D(j,j)=0;
    Temp_P(j,j)=0;
end

SDA_detail_import(1,:)=sum(Temp_E,1);
SDA_detail_import(2,:)=sum(Temp_L,1);
SDA_detail_import(3,:)=sum(Temp_S,1);
SDA_detail_import(4,:)=sum(Temp_D,1);
SDA_detail_import(5,:)=sum(Temp_P,1);

Result_driving{(i-1)*n_driving*49+1:(i-1)*n_driving*49+n_driving*49,6}=SDA_detail_import(:);


SDA_detail_export(1,:)=sum(Temp_E,2);
SDA_detail_export(2,:)=sum(Temp_L,2);
SDA_detail_export(3,:)=sum(Temp_S,2);
SDA_detail_export(4,:)=sum(Temp_D,2);
SDA_detail_export(5,:)=sum(Temp_P,2);

Result_driving{(i-1)*n_driving*49+1:(i-1)*n_driving*49+n_driving*49,7}=SDA_detail_export(:);


end


writetable(Result_driving,'Result_product.xlsx','Sheet',2);
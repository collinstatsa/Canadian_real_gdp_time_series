%% Matlab coding for canadian gdp analysis
clear all % efface toutes les variables en mémoire

diary 2022.out % sauvegarde résultats dans un fichier

% read the data

starts= xlsread('devoir1.xlsx','c1:c243');

dates = datetime('jan-1961')+calquarters(0:length(starts)-1);

%% a) Le graphique de la série

lstarts = log(starts);
diffstarts=diff(lstarts);
dattes = dates(2:end);
plot(dattes,diffstarts);
title({'Taux croissance pib réel canadien - 1961-2021','Désaisonnalisées au taux annuel'});
hold on;
yline(0)
hold off;

%% b) Les 24 premières autocorrélations et autocorrélations partielles
figure;autocorr(diffstarts,24);
figure;parcorr(diffstarts,24);

%% c) AIC et BIC pour le choix du modèle ARMA

% Estimation d'un modèle ARMA(1,1)
model1 = arima(1,0,1);
estmodel1=estimate(model1,diffstarts);
[res1,~,logL] = infer(estmodel1,diffstarts);

% Figure des résidus
model2 = arima(1,0,1);
estmodel2=estimate(model2,lstarts);

[res2,~,logL] = infer(estmodel2,lstarts);

plot(res2)

figure;
subplot(2,1,1);autocorr(res2,20)
subplot(2,1,2);parcorr(res2,20)


% Critère de minimisation

pMax = 3;
qMax = 3;

LogL = zeros(pMax+1,qMax+1);
SumPQ = LogL;

for p = 0:pMax
    for q = 0:qMax
        [p,q]
        Mdl = arima(p,0,q);
        [~,~,LogL(p+1,q+1)] = estimate(Mdl,diffstarts,'Display','off');
        SumPQ(p+1,q+1) = p+q;
    end
end

logL = reshape(LogL,(pMax+1)*(qMax+1),1); % convertit matrice en vecteur
numParams = reshape(SumPQ,(pMax+1)*(qMax+1),1)+1 ;
[aic,bic] = aicbic(logL,numParams,length(diffstarts));
AIC = reshape(aic,pMax+1,qMax+1)
minAIC = min(aic);
[bestP,bestQ] = find(AIC == minAIC);

disp('Modele choisi par AIC'); [bestP-1,bestQ-1]

BIC = reshape(bic,pMax+1,qMax+1)
minBIC = min(bic);
 [bestPbic,bestQbic]= find(BIC == minBIC);

disp('Modele choisi par BIC'); [bestPbic-1,bestQbic-1]

%% d) Reprise des sous-questions b) et c) pour la période 1961:2 à 2019:4

% Reshape de la base de données
dattes3=dates(2:236)
diffstarts3=diffstarts(2:236)

figure;autocorr(diffstarts3,24);
figure;parcorr(diffstarts3,24);

% Estimation d'un modèle ARMA(1,2)
model1 = arima(1,0,2);
estmodel1=estimate(model1,diffstarts3);
[res1,~,logL] = infer(estmodel1,diffstarts3);

%Autocorrélogramme des résidus du ARMA(1,2)

model2 = arima(1,0,2);
estmodel2=estimate(model2,lstarts);

[res2,~,logL] = infer(estmodel2,lstarts);

plot(res2)

figure;
subplot(2,1,1);autocorr(res2,20)
subplot(2,1,2);parcorr(res2,20)

% Critère de minimisation

pMax = 3;
qMax = 3;

LogL = zeros(pMax+1,qMax+1);
SumPQ = LogL;

for p = 0:pMax
    for q = 0:qMax
        [p,q]
        Mdl = arima(p,0,q);
        [~,~,LogL(p+1,q+1)] = estimate(Mdl,diffstarts3,'Display','off');
        SumPQ(p+1,q+1) = p+q;
    end
end

logL = reshape(LogL,(pMax+1)*(qMax+1),1); % convertit matrice en vecteur
numParams = reshape(SumPQ,(pMax+1)*(qMax+1),1)+1 ;
[aic,bic] = aicbic(logL,numParams,length(diffstarts3));
AIC = reshape(aic,pMax+1,qMax+1)
minAIC = min(aic);
[bestP,bestQ] = find(AIC == minAIC);

disp('Modele choisi par AIC'); [bestP-1,bestQ-1]

BIC = reshape(bic,pMax+1,qMax+1)
minBIC = min(bic);
 [bestPbic,bestQbic]= find(BIC == minBIC);

disp('Modele choisi par BIC'); [bestPbic-1,bestQbic-1]


%% e) Redefinir le taux de croissance puis refaire b) et c)

%%pibcrt=diff(starts)
Y_star = diff(starts)./starts(1:end-1)  %Calcul de Y_star

figure;autocorr(Y_star,24);
figure;parcorr(Y_star,24);
%%%Je dois modifier mes valeurs ici pour tenir compte de Y_star
% Estimation d'un modèle ARMA(1,1)
model1 = arima(1,0,1);
estmodel1=estimate(model1,Y_star);
[res1,~,logL] = infer(estmodel1,Y_star);

% Critère de minimisation

pMax = 3;
qMax = 3;

LogL = zeros(pMax+1,qMax+1);
SumPQ = LogL;

for p = 0:pMax
    for q = 0:qMax
        [p,q]
        Mdl = arima(p,0,q);
        [~,~,LogL(p+1,q+1)] = estimate(Mdl,Y_star,'Display','off');
        SumPQ(p+1,q+1) = p+q;
    end
end

logL = reshape(LogL,(pMax+1)*(qMax+1),1); % convertit matrice en vecteur
numParams = reshape(SumPQ,(pMax+1)*(qMax+1),1)+1 ;
[aic,bic] = aicbic(logL,numParams,length(Y_star));
AIC = reshape(aic,pMax+1,qMax+1)
minAIC = min(aic);
[bestP,bestQ] = find(AIC == minAIC);

disp('Modele choisi par AIC'); [bestP-1,bestQ-1]

BIC = reshape(bic,pMax+1,qMax+1)
minBIC = min(bic);
 [bestPbic,bestQbic]= find(BIC == minBIC);

disp('Modele choisi par BIC'); [bestPbic-1,bestQbic-1]








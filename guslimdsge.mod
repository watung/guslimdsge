
/*
 * Konvensi penamaan variabel adalah sebagai berikut:
 * - indeks yang dituliskan sebagai superscript dituliskan dengan
 *   underscore, misalnya: C_P, r_D, dan sebagainya.
 * - indeks yang dituliskan sebagai subscript dituliskan tanpa
 *   underscore, misalnya: epsilonu, PR, PK, dan sebagainya.
 */

var				epsilonu		//Intertemporal shock 
				epsilonh
				epsilonr
				epsilonrho
				epsilonR		//Housing preference
				epsilonL		//Labour preference
				epsilon_D
				epsilon_IR
				epsilon_E
				epsilon_IC
				psiut
				r
				rstar			//Suku bunga luar negeri
				taustar			//Pajak atas aliran modal masuk  
				r_D			//Suku bunga deposito
				r_K			//Suku bunga kapital 
				P			//Harga final 
				Pstar
				Bstar			//Pinjaman luar negeri
				Btotstar
				BGstar
				A
				C
				D			//Total deposito 
				E
				UCaksen
				PR			//Harga rumah
				PK			//Harga barang kapital 
				pi			//Inflasi
				piR
				pistar			//Inflasi luar negeri
				rho			//Premi risiko
				s			//Nilai tukar riil 
				e			//Nilai tukar nominal 
				y
				y_E
				ystar
				ytilde
				eta
				etastar
				mu
				Gamma
				B
				BPstar
				BEstar
				G
				B_E
				upsilon
				upsilonB
				lambdaP
				R_D
				R_B
				


				/* Variabel spesifik PATIENT HOUSEHOLD */			
				C_P			//Tingkat konsumsi pada Patient Household
				R_P			//Kepemilikan asset perumahan pada Patient Household
				T_P			//Pajak yang dikeluarkan Patient Household 						
				K_P			//Kapital pada Patient Household
				BG_P			//Obligasi pemerintah pada Patient Household
				Pi_P			//Total dividen pada Patient Household
				W_P			//Pendapatan dari tenaga kerja pada Patient Household
				L_P			//Jumlah jam bekerja pada Patient Household
				XP
				YP
				piWP

				/* Variabel spesifik IMPATIENT HOUSEHOLD */			
				C_I			//Tingkat konsumsi pada Impatient Household
				R_I			//Kepemilikan asset perumahan pada Impatient Household
				T_I			//Pajak yang dikeluarkan Impatient Household 
				W_I			//Pendapatan dari tenaga kerja pada Impatient Household
				L_I			//Jumlah jam bekerja pada Impatient Household
				B_IC
				B_IR
				m_I
				r_IR
				r_IC
				b_IC
				piWI
				XI
				YI
				
				/* Variabel spesifik ENTREPENEURS */			
				C_E			//Tingkat konsumsi pada Entrepeneurs
				P_E
				K_E
				r_E
				m_E
				
				/* Variabel spesifik DOMESTIC RETAILERS */			
				yH
				piH
				XH
				YH
				PH
				Pi_H

				/* Variabel spesifik IMPORTIR RETAILERS */			
				yF			
				piF
				XF
				YF
				PF
				PFstar
				Pi_F
				
				/* Variabel spesifik EXPORTING RETAILERS */			
				yHstar
				piHstar
				XHstar
				YHstar
				PHstar
				Pi_Hstar
				
				/* Variabel spesifik CAPITAL GOODS PRODUCERS */
				K
				epsiloniK
				iK

				/* Variabel spesifik HOUSING GOODS PRODUCERS */
				R
				epsiloniR
				iR
				
				/* Variabel spesifik BANK */
				K_B
				Pi_B
				BG_B;
				

			
		
parameters			xi			//Tingkat external habit formation
				sigmaC			//Elastisitas intertemporal dari substitusi di konsumsi
				sigmaR			//Elastisitas intertemporal dari substitusi di rumah
				sigmaL			//Elastisitas dari bekerja
				sigmaE
				alphaTD
				alphaTW
				alphaTG
				alphaTK
				alphaTPi
				alphaE
				deltaR			//Tingkat depresiasi perumahan
				deltaK			//Tingkat depresiasi kapital
				deltaE
				deltaB
				w_B
				epsilonw
				epsilonW
				kappaK
				kappaKB
				kappaR
				kappaD
				kappaE
				kappaIR
				kappaIC
				omegaB
				rbar
				pibar
				PhiR
				Phipi
				Phiy
				nu

				/* Parameter spesifik pada PATIENT HOUSEHOLD */	
				betaP  			//Discount factor pada Patient Household
				thetaWP			
				thetaP
				gamma_P
				
				/* Parameter spesifik pada IMPATIENT HOUSEHOLD */	
				betaI			//Discount factor pada Impatient Household
				thetaWI
				thetaI
				gamma_I
				
				/* Parameter spesifik pada ENTREPENEURS */	
				betaE			//Discount factor pada Entrepreneurs
				gamma_E

				/* Parameter spesifik pada DOMESTIC RETAILERS */			
				epsilonH
				thetaH
				
				/* Parameter spesifik pada IMPORTIR RETAILERS */			
				epsilonF
				thetaF

				/* Parameter spesifik pada EXPORTING RETAILERS */			
				epsilonHstar
				thetaHstar
				muHstar;


model;


/*----------------------------*/
/*    I. PATIENT HOUSEHOLD    */
/*----------------------------*/

/* Dari persamaan PI.1 */
(betaP*epsilonu(+1)*(C_P(+1)-xi*C_P)^(-sigmaC)/P(+1))*((1+r_D*(1-alphaTD))/pi(+1))=epsilonu*(C_P-xi*C_P(-1))^(-sigmaC)/P;

/* Dari persamaan PI.2 */
epsilonu*epsilonR*R_P^(-sigmaR)+betaP*epsilonu(+1)*(C_P(+1)-xi*C_P)^(-sigmaC)/P(+1)*PR(+1)*(1-deltaR)=epsilonu*(C_P-xi*C_P(-1))^(-sigmaC)/P*PR;

/* Dari persamaan PI.3 */
D=W_P*L_P+(1+r_D(-1))*D(-1)/pi+(1+r(-1))*BG_P(-1)/pi+r_K*PK*K_P(-1)+e*Bstar+Pi_P*P*C_P-PK*(K_P-(1-deltaK)*K_P(-1))-PR*(R_P-(1-deltaR)*R_P(-1))-BG_P-e*(1+rho(-1))*(1+rstar(-1))*(1+taustar(-1))*Bstar(-1)/pi-(alphaTW*W_P*L_P+alphaTD*r_D(-1)*D(-1)+alphaTG*r(-1)*BG_P(-1)+alphaTK*r_K*PK*K_P(-1)+alphaTPi*Pi_P);

/* Dari persamaan PI.4 */
Pi_P=Pi_H+Pi_Hstar+Pi_F+(1-w_B)*Pi_B;

/* Dari persamaan PI.5 */
T_P=alphaTW*W_P*L_P+alphaTD*r_D(-1)*D(-1)+alphaTG*r(-1)*BG_P(-1)+alphaTK*r_K*PK*K_P(-1)+alphaTPi*Pi_P;

/* Dari persamaan PI.6 */
(1+r_D*(1-alphaTD))/((1+rho)*(1+rstar)*(1+taustar))*E*(pistar(+1)/pi(+1))=E*(s(+1)/s);

/* Dari persamaan PI.7 */
pi=P/P(-1);

/* Dari persamaan PI.8 */
s=e*Pstar/P;

/* Dari persamaan PI.9 */
XP=epsilonu*(piWP^(epsilonw-1)*W_P*L_P/P*UCaksen)+betaP*thetaWP*XP(+1);

/* Dari persamaan PI.10 */
YP=epsilonu*(epsilonL*(piWP^epsilonw*L_P)^(1+sigmaL))+betaP*thetaWP*YP(+1);

/* Dari persamaan PI.11 */
(1-alphaTW)*(epsilonw-1)*((piWP^(1-epsilonW)-thetaP*piWP(-1)^(1-epsilonW))/(1-thetaP))^((epsilonw*sigmaL+1)/(1-epsilonW))/(epsilonw*gamma_P^(-sigmaL))*XP=YP;


/*-------------------------------*/
/*    II. IMPATIENT HOUSEHOLD    */
/*-------------------------------*/

/* Dari persamaan PII.1 */
C_I=W_I*L_I*(1-alphaTW)/P+B_IC/P+m_I*E*(PR(+1)*(1-deltaR)*R_I)/(P*(1+r_IR))-(1+r_IC(-1))/P*b_IC(-1)/pi-m_I(-1)*E(-1)*(PR*(1-deltaR)*R_I(-1))/(P*pi)-PR*(R_I-(1-deltaR)*R_I(-1))/P;

/* Dari persamaan PII.2 */
epsilonu*epsilonh*R_I^(-sigmaR)+epsilonu*(C_I-xi*C_I(-1))^(-sigmaC)*(m_I*E*(PR(+1)*(1-deltaR)/(P*(1+r_IR)))-PR/P)+betaI*epsilonu*(C_I(+1)-xi*C_I)^(-sigmaC)*(betaI*epsilonu*(C_I(+1)-xi*C_I)^(-sigmaC)*PR(+1)*(1-deltaR)/P(+1)-betaI*epsilonu*m_I*E*PR(+1)*(1-deltaR)/(P(+1)*pi(+1)))=0;

/* Dari persamaan PII.3 */
B_IR=m_I*E*(PR(+1)*(1-deltaR)*R_I)/(1+r_IR);

/* Dari persamaan PII.4 */
T_I=alphaTW*W_I*L_I;

/* Dari persamaan PII.5 */
betaI*epsilonu(+1)*(C_I(+1)-xi*C_I)^(-sigmaC)/P(+1)*(1+r_IC)/pi(+1)=epsilonu*(C_I-xi*C_I(-1))^(-sigmaC)/P;

/* Dari persamaan PII.6 */
XI=epsilonu*(piWI^(epsilonw-1)*W_I*L_I/P*UCaksen)+betaI*thetaWI*XI(+1);

/* Dari persamaan PII.7 */
YI=epsilonu*(epsilonL*(piWI^epsilonw*L_I)^(1+sigmaL))+betaI*thetaWI*YI(+1);

/* Dari persamaan PII.8 */
(1-alphaTW)*(epsilonw-1)*((piWI^(1-epsilonW)-thetaI*piWI(-1)^(1-epsilonW))/(1-thetaI))^((epsilonw*sigmaL+1)/(1-epsilonW))/(epsilonw*gamma_I^(-sigmaL))*XI=YI;


/*-------------------------*/
/*    III. ENTREPENEURS    */
/*-------------------------*/

/* Dari persamaan PIII.1 */
betaE*epsilonu(+1)*(C_E(+1)-xi*C_E)^(-sigmaC)/P(+1)*PK(+1)*(1-deltaE)+betaE*epsilonu(+1)*(C_E(+1)-xi*C_E)^(-sigmaC)/P(+1)*P_E(+1)*A(+1)*sigmaE*alphaE*K_E^(sigmaE*alphaE-1)*(K_P^(1-sigmaE))^alphaE*(L_P(+1)^gamma_E*L_I(+1)^(1-gamma_E))^(1-alphaE)+epsilonu*(C_E-xi*C_E(-1))^(-sigmaC)/P*m_E*E*PK(+1)*(1-deltaK)/(1+r_E)=epsilonu*(C_E-xi*C_E(-1))^(-sigmaC)/P*PK+betaE*epsilonu(+1)*(C_E(+1)-xi*C_E)^(-sigmaC)/P(+1)*m_E*E*PK(+1)*(1-deltaK)/pi(+1);

/* Dari persamaan PIII.2 */
r_K(+1)*PK(+1)=P_E(+1)*A(+1)*(1-sigmaE)*alphaE*K_E^(sigmaE*alphaE)*K_P^((1-sigmaE)*alphaE-1)*(L_P(+1)^gamma_E*L_I(+1)^(1-gamma_E))^(1-alphaE);

/* Dari persamaan PIII.3 */
betaE*epsilonu(+1)*(C_E(+1)-xi*C_E)^(-sigmaC)/P(+1)*e(+1)*(1+rho)*(1+rstar)*(1+taustar)/pi(+1)=epsilonu*(C_E-xi*C_E(-1))^(-sigmaC)/P;


/*------------------------------*/
/*    IV. DOMESTIC RETAILERS    */
/*------------------------------*/

/* Dari persamaan PIV.1 */
XH=yH*piH^(epsilonH-1)+betaP*thetaH*XH(+1);

/* Dari persamaan PIV.2 */
YH=yH*piH^epsilonH*P_E/PH+betaP*thetaH*YH(+1);

/* Dari persamaan PIV.3 */
(1/(1-thetaH)*piH^(1-epsilonH)-thetaH*piH(-1)^(1-epsilonH))^(1/(1-epsilonH))*E*XH=epsilonH/(epsilonH-1)*E*YH;

/* Dari persamaan PIV.4 */
piH=PH/PH(-1);

/* Dari persamaan PIV.5 */
yH=eta*(PH/P)^(-((1+mu)/mu))*y;

/* Dari persamaan PIV.6 */
Pi_H=(PH-P_E)*yH;


/*------------------------------*/
/*    IV. IMPORTIR RETAILERS    */
/*------------------------------*/

/* Dari persamaan PV.1 */
XF=yF*piF^(epsilonF-1)+betaP*thetaF*XF(+1);

/* Dari persamaan PV.2 */
YF=yF*piF^epsilonF*e*PFstar/PF+betaP*thetaF*YF(+1);

/* Dari persamaan PV.3 */
(1/(1-thetaF)*piF^(1-epsilonF)-thetaF*piF(-1)^(1-epsilonF))^(1/(1-epsilonF))*E*XF=epsilonF/(epsilonF-1)*E*YF;

/* Dari persamaan PIV.4 */
piF=PF/PF(-1);

/* Dari persamaan PIV.5 */
yH=(1-eta)*(PF/P)^(-((1+mu)/mu))*y;

/* Dari persamaan PIV.6 */
Pi_F=(PF-e*PFstar)*yF;


/*-------------------------------*/
/*    VI. EXPORTING RETAILERS    */
/*-------------------------------*/

/* Dari persamaan PVI.1 */
XHstar=yHstar*piHstar^(epsilonHstar-1)+betaP*thetaHstar*XHstar(+1);

/* Dari persamaan PVI.2 */
YHstar=yHstar*piHstar^epsilonHstar*P_E/PHstar+betaP*thetaHstar*YHstar(+1);

/* Dari persamaan PVI.3 */
(1/(1-thetaHstar)*piHstar^(1-epsilonHstar)-thetaHstar*piHstar(-1)^(1-epsilonHstar))^(1/(1-epsilonHstar))*E*XHstar=epsilonHstar/(epsilonHstar-1)*E*YHstar;

/* Dari persamaan PVI.4 */
piHstar=PHstar/PHstar(-1);

/* Dari persamaan PVI.5 */
yHstar=(1-etastar)*(PHstar/Pstar)^(-((1+muHstar)/muHstar))*ystar;

/* Dari persamaan PVI.6 */
Pi_Hstar=(e*PHstar-P_E)*yHstar;


/*----------------------------------*/
/*    VII. FINAL GOODS PRODUCERS    */
/*----------------------------------*/

/* Dari persamaan PVII.1 */
P^(-1/mu)=eta*PH^(-1/mu)+(1-eta)*PF^(-1/mu);

/* Dari persamaan PVII.2 */
yH=eta*(PH/P)^(-(1+mu)/mu)*y;

/* Dari persamaan PVII.3 */
yF=eta*(PF/P)^(-(1+mu)/mu)*y;


/*-------------------------------------*/
/*    VIII. CAPITAL GOODS PRODUCERS    */
/*-------------------------------------*/

/* Dari persamaan PVIII.1 */
1/kappaK+betaP*PK(+1)/PK*epsiloniK(+1)/epsiloniK*(iK(+1)/iK)^2*(iK(+1)/iK-1)=(1/2*(iK/iK(-1)-1)^2+iK/iK(-1)*(iK/iK(-1)-1))+P/(kappaK*PK*epsiloniK);

/* Dari persamaan PVIII.2 */
K=(1-deltaK)*K(-1)+epsiloniK*(1-1/2*kappaK*(iK/iK(-1)-1)^2)*iK;


/*-----------------------------------*/
/*    IX. HOUSING GOODS PRODUCERS    */
/*-----------------------------------*/

/* Dari persamaan PIX.1 */
1/kappaR+betaP*PR(+1)/PR*epsiloniR(+1)/epsiloniR*(iR(+1)/iR)^2*(iR(+1)/iR-1)=(1/2*(iR/iR(-1)-1)^2+iR/iR(-1)*(iR/iR(-1)-1))+P/(kappaR*PR*epsiloniR);

/* Dari persamaan PIX.2 */
R=(1-deltaR)*R(-1)+epsiloniR*(1-1/2*kappaR*(iR/iR(-1)-1)^2)*iR;

/* Dari persamaan PIX.3 */
piR=PR/PR(-1);


/*-----------------*/
/*     X. BANK     */
/*-----------------*/

/* Dari persamaan PX.1 */
K_B=(1-deltaB)*K_B(-1)+w_B*Pi_B(-1);

/* Dari persamaan PX.2 */
BG_B=(1-Gamma)*D+K_B-B;

/* Dari persamaan PX.3 */
B=B_E+B_IC+B_IR;

/* Dari persamaan PX.4 */
Pi_B=r_E*B_E+r_IR*B_IR+r_IC*B_IC+r*BG_B-r_D*D-kappaKB/2*(B/(omegaB*K_B)-upsilon)^2*K_B-kappaD/2*(r_D/r_D(-1)-1)^2*r_D*D-kappaE/2*(r_E/r_E(-1)-1)^2*r_E*B_E-kappaIR/2*(r_IR/r_IR(-1)-1)^2*r_IR*B_IR-kappaIC/2*(r_IC/r_IC(-1)-1)^2*r_IC*B_IC;

/* Dari persamaan PX.5 */
r*(1-Gamma)=R_D;

/* Dari persamaan PX.6 */
R_B=r-omegaB*kappaKB*(K_B/(omegaB*B)-upsilonB)*(K_B/(omegaB*B))^2;

/* Dari persamaan PX.7 */
epsilon_D+betaP*lambdaP(+1)/lambdaP*kappaD*(r_D(+1)/r_D-1)*r_D(+1)^2/r_D^2*D(+1)/D=1+epsilon_D*R_D/r_D+kappaD*(r_D/r_D(-1)-1)*r_D/r_D(-1);

/* Dari persamaan PX.8 */
1+epsilon_IR*R_B/r_IR+betaP*lambdaP(+1)/lambdaP*kappaIR*(r_IR(+1)/r_IR-1)*r_IR(+1)^2/r_IR^2*B_IR(+1)/B_IR=epsilon_IR+kappaIR*(r_IR/r_IR(-1)-1)*r_IR/r_IR(-1);

/* Dari persamaan PX.9 */
1+epsilon_IC*R_B/r_IC+betaP*lambdaP(+1)/lambdaP*kappaIC*(r_IC(+1)/r_IC-1)*r_IC(+1)^2/r_IC^2*B_IC(+1)/B_IC=epsilon_IC+kappaIC*(r_IC/r_IC(-1)-1)*r_IC/r_IC(-1);

/* Dari persamaan PX.10 */
1+epsilon_E*R_B/r_E+betaP*lambdaP(+1)/lambdaP*kappaE*(r_E(+1)/r_E-1)*r_E(+1)^2/r_E^2*B_E(+1)/B_E=epsilon_E+kappaE*(r_E/r_E(-1)-1)*r_E/r_E(-1);


/*------------------------*/
/*    XI. CENTRAL BANK    */
/*------------------------*/

/* Dari persamaan PXI.1 */
(1+r)=(1+rbar)^(1-PhiR)*(1+r(-1))^PhiR*((pi/pibar)^Phipi*(y/y(-1))^Phiy)^(1-PhiR)*epsilonr;


/*-----------------------*/
/*    XII. GOVERNMENT    */
/*-----------------------*/

/* Dari persamaan PXII.1 */
P*G+e*(1+rho(-1))*(1+rstar(-1))*(1+taustar(-1))*BGstar(-1)+(1+r(-1))*BG_P(-1)+(1+r(-1))*BG_B(-1)=(T_P+T_I+e*taustar(-1)*BGstar)+e*BGstar+BG_P+BG_B;


/*----------------------------------------*/
/*    XIII. MARKET CLEARING CONDITIONS    */
/*----------------------------------------*/

/* Dari persamaan PXIII.1 */
y_E=yH+yHstar;

/* Dari persamaan PXIII.2 */
y=C+iK+iR+G+psiut*K(-1);

/* Dari persamaan PXIII.3 */
C=gamma_P*C_P+gamma_I*C_I+gamma_E*C_E;

/* Dari persamaan PXIII.4 */
R=gamma_P*R_P+gamma_I*R_I;

/* Dari persamaan PXIII.5 */
K=gamma_P*K_P+gamma_E*K_E;

/* Dari persamaan PXIII.6 */
P*ytilde=P*y+e*PHstar*yHstar-e*PFstar*yF;

/* Dari persamaan PXIII.7 */
e*PFstar*yF+e*(1+rho(-1))*(1+rstar(-1))*Btotstar(-1)=e*PHstar*yHstar+e*Btotstar;

/* Dari persamaan PXIII.8 */
Btotstar=BPstar+BEstar+BGstar;

/* Dari persamaan PXIII.9 */
(1+rho)=exp(-nu*(e*Btotstar/(P*ytilde)))*epsilonrho;



end;

function     markov(res,Vs,B)

if nargin == 0
%Definisco res vettore delle resistenze arbitrario se non dato in ingresso
    res=[10 5 2 1 2 10 5 10];
end

if nargin < 2
%Definisco Vs arbitrario se non dato in ingresso.
    Vs=1;
end
if nargin < 3
%Se non viene data una matrice per legge di kirchoff assegno quella dell'esercizio 2.4 del Moler.
     B = [ 1 -1  0  0  0  
       1  0 -1  0  0 
       1  0  0 -1  0 
       0 -1  1  0  0 
       0  0  1 -1  0 
       0 -1  0  0  1 
       0  0 -1  0  1 
       0  0  0 -1  1 ]; 
end

% Qualche condizione per assicurarci che la matrice sia una matrice per
% legge delle correnti di Kirchoff:

l=length(nonzeros(B));
if sum(B')~=0 | abs(nonzeros(B))~= ones(l,1)
    error('La matrice inserita non può essere una matrice per legge di kirchoff delle correnti')
end


% % Definiamo la matrice B che descrive la legge di Kirchoff per le correnti
% % che coincide con quella dell'esercizio 2.4 del Moler. Le righe di tale
% % matrice rappresentano i collegamenti tra i vertici su cui è presente una 
% % resistenza mentre le colonne rappresentano i vertici

%Definiamo un vettore contenente i valori delle conduttanze. Bisogna
%prestare attenzione all'ordine degli elementi del vettore che andranno in 
% accordo con la matrice B ovvero il primo elemento di g sarà la 
% conduttanza tra i vertici 1 e 2, il secondo sarà la conduttanza 
% tra 1 e 3 etc...

    g=1./res; 

% Defianiamo la matrice delle conduttanze:

    G = B'*diag(g)*B;

% Attraverso un ciclo creiamo la matrice del processo markoviano che
% descrive il nostro circuito. Per definizione W(i,j) = I(i,j)-(G(i,j)/G(i,i))
% dove I rappresenta la matrice identità:

    I = ones(length(G),1);
    I=diag(I);
    i=1;
    j=1;
    while i < length(G)+1
        j=1;
        while j < length(G)+1
            W0(i,j) = I(i,j)-(G(i,j)/G(i,i));
            j=j+1;
        end
        i=i+1;
    end
    
    disp('matrice stocastica irriducibile')
    disp(W0)

% %la matrice stocastica così creata è bidirezionale per ogni ramo, ovvero è
% una matrice stocastica irriducibile. 
% Per ogni coppia di stati i,j posso andare da i a j e poi tornare indietro.
% Rappresenta in modo corretto la situazione del circuito senza generatore.
%L'inserimento del generatore fa si che si preferisca almeno una direzione
%ovvero quella che connette i vertici 3 e 5. E' necessario rendere tale
%nodo diretto. Scegliamo che un random walker non possa andare da 3 a 5 ma
%solo da 5 a 3. Questo viene fatto ponendo a 0 l'elemento di matrice
%corrispondente e normalizzando (si ridistribuiscono equamente le
%probabilità):

    W0(3,5)=0; %impedisco 3 -> 5
    W0(3,:)=W0(3,:)./sum(W0(3,:));%normalizzazione
    
% E' inoltre necessario modificare la matrice delle conduttanze tenendo
% conto che il nodo 3 non risente della conduttanza data da r35:

    G(3,3)=(1/res(2))+(1/res(4))+(1/res(5));
    disp('matrice stocastica con correzioni date dalla presenza del generatore')
    disp(W0)
    disp('isstochastic') %verifichiamo che la matrice è stocastica per righe
    disp(sum(W0'))
    pause
    g=digraph(W0); %creiamo il grafico diretto della matrice W0
    plot(g);
    
% Sappiamo che se W descrive un processo markoviano allora W.^n ci dice
% dove possiamo trovarci dopo n 'passi'. Supponendo n -> infinito (storia
% infinita) otterremo la soluzione di equilibrio statistico. Un teorema 
% importante in questo contesto è il teorema di Perron-Frobenius il quale
% afferma che esiste un autovalore massimo reale positivo con molteplicità 
% 1 e corrispondente autospazio di dimensione uno. Sapendo che
% le matrici stocastiche hanno spettro degli autovalori contenuto nel 
% cerchio unitario, |h|<=1 allora se la soluzione di equilibrio statistico
% peq (t.c. W*peq=peq) esiste essa è unica. Il teorema tuttavia si applica
% a matrici irriducibili e la presenza del generatore nel processo non
% soddisfa tale proprietà. Non è garantita la presenza della soluzione di
% equilibrio statistico.

    H=W0^10000;%Definiamo H come la matrice di un processo markoviano di 1e4 passi
    H=H';
    
%la soluzione di equilibrio esiste infatti le colonne della matrice
%W0^10000 sono tutte uguali. Basterà prenderne una.

    u0=H(:,1);%prendo la prima colonna di H
    u0(:,1)=u0(:,1)./u0(3,1);% la normalizzo dividendo per la terza entrata
    disp('Matrice stocastica dopo 10000 passi')
    disp(H)

%Calcoliamo invece esplicitamente l'autovettore di W0 corrispondente ad
%autovalore 1:
    [u,e]=eig(full(W0'));
    e=diag(abs(e)); %gli autovalori stanno sulla diagonale della matrice 'e' 
    disp('autovalori della matrice stocastica con correzione del generatore')
    disp(e)%è presente 1 tra gli autovalori 
    pause
    
 %l'autovettore corrispondente all'equilibrio è in prima posizione. 
% Lo normalizziamo dividendolo per la sua componente più grande.

    u(:,1)=u(:,1)./u(3,1);
    
%Utilizzo markovchain per simulare un random walk sul mio circuito. In f
%otteniamo le frequenze relative alle visite sui vertici del grafo.

    [~,f]=markovchain(g,1e5,1,-1);%eseguo 100000 passi partendo da 1 e con una pjump nulla
    f=f./f(3);%normalizzo ancora dividendo per la terza componente
    disp('Contronto i risultati ottenuti tramite rispettivamente markovchain.m, autovettore associato all autovalore 1 e colonna di W0^10000')
    disp([f u(:,1) u0])
    disp('Controlliamo la precisione tra le misure')
    disp('norma tra random walk e autovettore')
    disp(norm(f-u(:,1)))
    disp('norma tra random walk e soluzione di equilibrio dopo 10000 passi')
    disp(norm(f-u0))
    disp('norma tra soluzione di equilibrio dopo 10000 passi e autovettore')
    disp(norm(u(:,1)-u0))
    pause
    
% Trovo i potenziali sapendo che u(k), componente k-esima della 
% soluzione di equilibrio soddisfa u(k)=v(k)*G(k,k) dove v(k) è il
% potenziale sul vertice k.

    volt=diag(diag(u0)./G);%utilizzo l'operatore / di matlab e prendo i valori sulla diagonale

% Normalizzo il vettore dei potenziali in modo che il vertice 3 sia a
% potenziale 1. In tale modo la moltiplicazione del vettore per il
% voltaggio del generatore ci darà i valori corretti dei potenziali.
    
    volt=Vs*(volt./volt(3));
    disp('potenziali')
    disp(volt)

% Bisogna modificare la matrice B in modo che possa dare le differenze
% di potenziali corrette tra i vertici. Infatti, se tenessimo la B di
% partenza, avremmo una differenza V(5)-V(3) mentre la corrente che scorre 
% sulla resistenza r35 è semplicemente (V(5)-0)/r35.Di conseguenza la 
% nuova matrice non avrà la componente (7,3):
    
    B(7,3)=0;
    deltav=B*volt; %vettore contenente le differenze di potenziale
    I=diag(deltav./res);%vettore contenente le correnti
    disp('vettore contente le correnti sui rami nell ordine definito da B')
    disp(I)
    pause
    disp('Verifichiamo la legge di Kirchoff, la somma algebrica delle correnti su un vertice deve essere nulla')
    disp('vertice 1')
    disp(I(3)+I(2)+I(1))
    disp('vertice 2')
    disp(I(6)+I(4)+I(1))
    disp('vertice 3')
    disp(-I(4)-I(5)+I(2)+I(7))%La corrente uscente de essere uguale a quella che passa sulla resistenza r35
     disp('vertice 4')
    disp(I(8)+I(5)+I(3))
    disp('vertice 5')
    disp(I(8)+I(6)+I(7))
    
end
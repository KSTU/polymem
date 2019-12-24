%
pkg load io
arg_list=argv();
if(length(arg_list)<3)
	fprintf("Error. Usege: octave polymem.m [monomer file] [monomer chain] [numers of molecules] [x-y width] \n");
	break
else
mono_file_name=arg_list{1}
mono_chain_number=str2num(arg_list{2})
mol_number=str2num(arg_list{3})
box_h=str2num(arg_list{4})

endif
#ввод исходный данных
file_id=fopen(mono_file_name,"r")
	tempstring=fscanf(file_id,"%s",1);
	mono_atom_number=fscanf(file_id,"%f",1);
	temp_string=fgetl(file_id);
	for ii=1:mono_atom_number
		temp_string=fgetl(file_id);
		mono_struct_name(ii,:)= strtrim(substr(temp_string,6,5));
		mono_atom_name(ii,:)= strtrim(substr(temp_string,11,5));
		mono_atom_x(ii)=str2num(substr(temp_string,21,8));
		mono_atom_y(ii)=str2num(substr(temp_string,29,8));
		mono_atom_z(ii)=str2num(substr(temp_string,37,8));
#		
	endfor
fclose(file_id)

#перенос первого автома в начало координат
temp_x=mono_atom_x(1);
temp_y=mono_atom_y(1);
temp_z=mono_atom_z(1);
for ii=1:mono_atom_number
	mono_atom_x(ii)=mono_atom_x(ii)-temp_x;
	mono_atom_y(ii)=mono_atom_y(ii)-temp_y;
	mono_atom_z(ii)=mono_atom_z(ii)-temp_z;
endfor

#перенос координат вдоль оси z
temp_x=mono_atom_x(mono_atom_number);
temp_y=mono_atom_y(mono_atom_number);
temp_z=mono_atom_z(mono_atom_number);

%поворот по оси z

cosa=0;	%	temp_x/sqrt(temp_x^2+temp_y^2);
sina=sqrt(1-cosa^2)
cosb=0; %temp_z/sqrt(temp_y^2+temp_z^2);
sinb=sqrt(1-cosb^2)
cosg=temp_y/sqrt(temp_y^2+temp_x^2);
sing=sqrt(1-cosg^2)

for ii=1:mono_atom_number
	temp_x=mono_atom_x(ii)*(cosa*cosg-sina*cosb*sing)-(cosa*sing+sina*cosb*cosg)*mono_atom_y(ii)+(sina*sinb)*mono_atom_z(ii);
	temp_y=(sina*cosg+cosa*cosb*cosg)*mono_atom_x(ii)+(-sina*sing+cosa*cosb*cosg)*mono_atom_y(ii)-(cosa*sinb)*mono_atom_z(ii);
	temp_z=sina*sing*mono_atom_x(ii)+sinb*cosg*mono_atom_y(ii)+cosb*mono_atom_z(ii);
	mono_atom_x(ii)=temp_x;
	mono_atom_y(ii)=temp_y;
	mono_atom_z(ii)=temp_z;
endfor

mono_atom_x
mono_atom_y
mono_atom_z

%поворот по оси x
temp_x=mono_atom_x(mono_atom_number);
temp_y=mono_atom_y(mono_atom_number);
temp_z=mono_atom_z(mono_atom_number);

cosa=0;	%	temp_x/sqrt(temp_x^2+temp_y^2);
sina=sqrt(1-cosa^2)
cosb=temp_z/sqrt(temp_x^2+temp_z^2);
sinb=sqrt(1-cosb^2)
cosg=0;	%temp_y/sqrt(temp_y^2+temp_x^2);
sing=sqrt(1-cosg^2)

for ii=1:mono_atom_number
	temp_x=mono_atom_x(ii)*(cosa*cosg-sina*cosb*sing)-(cosa*sing+sina*cosb*cosg)*mono_atom_y(ii)+(sina*sinb)*mono_atom_z(ii);
	temp_y=(sina*cosg+cosa*cosb*cosg)*mono_atom_x(ii)+(-sina*sing+cosa*cosb*cosg)*mono_atom_y(ii)-(cosa*sinb)*mono_atom_z(ii);
	temp_z=sina*sing*mono_atom_x(ii)+sinb*cosg*mono_atom_y(ii)+cosb*mono_atom_z(ii);
	mono_atom_x(ii)=temp_x;
	mono_atom_y(ii)=temp_y;
	mono_atom_z(ii)=temp_z;
endfor

mono_atom_x
mono_atom_y
mono_atom_z

bond_l=0.154;

% получение 
id=0;
for ii=1:mono_chain_number
	for jj=1:mono_atom_number
		id++;
		poly_atom_x(id)=mono_atom_x(jj);
		poly_atom_y(id)=mono_atom_y(jj);
		poly_atom_z(id)=mono_atom_z(jj)+(mono_atom_z(mono_atom_number)+bond_l)*(ii-1);
		poly_struct_name(id,:)=mono_struct_name(jj,:);
		poly_atom_name(id,:)=mono_atom_name(jj,:);
	endfor
endfor
poly_atom_number=mono_atom_number*mono_chain_number;

#poly gro out
file_id=fopen("poly.gro","w");
fprintf(file_id,"polymer\n");
fprintf(file_id,"%5d\n", columns(poly_atom_x))
for ii=1:columns(poly_atom_x)
	fprintf(file_id,"%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",1,poly_struct_name(ii,:),poly_atom_name(ii,:),ii,poly_atom_x(ii),poly_atom_y(ii),poly_atom_z(ii),0,0,0)
endfor
fprintf(file_id," %4.1f  %4.1f  %4.1f",0.5,0.5,max(poly_atom_z)+0.1)
fclose(file_id);

#read structure
file_id=fopen("top.in","r")
bond_n=fskipl(file_id,Inf)
frewind(file_id)
for ii=1:bond_n
	bond_atom1(ii)=fscanf(file_id,"%f",1);
	bond_atom2(ii)=fscanf(file_id,"%f",1);
endfor
fclose(file_id)

id=0
#creat poly bonds
for ii=1:mono_chain_number
	for jj=1:bond_n
		id++;
		poly_bond_1(id)=bond_atom1(jj)+(ii-1)*mono_atom_number;
		poly_bond_2(id)=bond_atom2(jj)+(ii-1)*mono_atom_number;
	endfor
endfor
for ii=1:mono_chain_number-1
	id++;
	poly_bond_1(id)=mono_atom_number*ii
	poly_bond_2(id)=mono_atom_number*ii+1
endfor
poly_bond_number=id
#get bond matrix
bond_matrix=zeros(poly_atom_number,poly_atom_number);
for ii=1:poly_bond_number
	bond_matrix(poly_bond_1(ii),poly_bond_2(ii))=1;
	bond_matrix(poly_bond_2(ii),poly_bond_1(ii))=1;
endfor
bond_matrix
#get angles
id_an=0;
poly_atom_number
for ii=1:poly_atom_number
	%ii
	for jj=ii+1:poly_atom_number
		%jj
		if (bond_matrix(ii,jj)==1)
			for kk=jj+1:poly_atom_number
				%kk
				if(bond_matrix(jj,kk)==1)
#					bond_matrix(ii,jj)
#					bond_matrix(jj,kk)
					id_an++;
					poly_angle_1(id_an)=ii;
					poly_angle_2(id_an)=jj;
					poly_angle_3(id_an)=kk;
				endif
			endfor
		endif
	endfor
endfor
poly_angle_number=id_an
#get dihedrals
id_dih=0;
for ii=1:poly_atom_number
	for jj=ii+1:poly_atom_number
		if (bond_matrix(ii,jj)==1)
			for kk=jj+1:poly_atom_number
				if (bond_matrix(jj,kk)==1)
					for ll=kk+1:poly_atom_number
						if (bond_matrix(kk,ll)==1)
							id_dih++;
							poly_dih1(id_dih)=ii;
							poly_dih2(id_dih)=jj;
							poly_dih3(id_dih)=kk;
							poly_dih4(id_dih)=ll;
						endif
					endfor
				endif
			endfor
		endif
	endfor
endfor
poly_dih_number=length(poly_dih1)
#write structure

file_id=fopen("out.top","w");
fprintf(file_id,"[ defaults ]\n")
fprintf(file_id,"1  2 no  1.0   1.0\n")
fprintf(file_id,"[ atomtypes ]\n")
fprintf(file_id,"C   14.0000   0.00 A    0.395   0.382444\n")
fprintf(file_id,"[ bondtypes ]\n")
fprintf(file_id,"C    C      1        0.154    500.00000\n")
fprintf(file_id,"[ angletypes ]\n")
fprintf(file_id,"C    C    C      1  114.0    530.2295918367\n")
fprintf(file_id,"[ dihedraltypes ]\n")
fprintf(file_id,"C    C    C    C      5  0.0  3.0119585918   -0.5785016939  6.7133004898   \n")
fprintf(file_id," [ moleculetype ]\n")
fprintf(file_id,"POL  3\n")
fprintf(file_id,"[ atoms ]\n")
for ii=1:poly_atom_number
	fprintf(file_id,"%d  C      1 Pol C    1   0.0  1.0\n",ii)
endfor

fprintf(file_id,"[ bonds ]\n")
for ii=1:poly_bond_number
	fprintf(file_id,"%d %d %d\n",poly_bond_1(ii),poly_bond_2(ii),1)
endfor


fprintf(file_id,"[ pairs ]\n")
fprintf(file_id,"[ angles ]\n")
for ii=1:poly_angle_number
	fprintf(file_id,"%d %d %d 1\n",poly_angle_1(ii),poly_angle_2(ii),poly_angle_3(ii))
endfor
fprintf(file_id,"[ dihedrals ]\n")
for ii=1:poly_dih_number
	fprintf(file_id,"%d %d %d %d 5\n",poly_dih1(ii),poly_dih2(ii),poly_dih3(ii),poly_dih4(ii))
endfor
fprintf(file_id,"\n[ system ]\n")
fprintf(file_id,"temp\n")
fprintf(file_id,"[ molecules ]\n")
fprintf(file_id,"POL %d",mol_number)

fclose(file_id);

#membrane gro out
file_id=fopen("membrane.gro","w");
fprintf(file_id,"membrane\n");
fprintf(file_id,"%5d\n", columns(poly_atom_x)*mol_number)
id=0;
chk=ceil(sqrt(mol_number))
for ii=1:chk
	for jj=1:chk
		id++;
		deltax(id)=(ii-0.5)*box_h/chk;
		deltay(id)=(jj-0.5)*box_h/chk;
	endfor
endfor

for jj=1:mol_number
	for ii=1:columns(poly_atom_x)
		fprintf(file_id,"%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",jj,poly_struct_name(ii,:),poly_atom_name(ii,:),ii+(jj-1)*columns(poly_atom_x),poly_atom_x(ii)+deltax(jj),poly_atom_y(ii)+deltay(jj),poly_atom_z(ii),0,0,0)
	endfor
endfor
fprintf(file_id," %4.1f  %4.1f  %4.1f ",box_h,box_h,max(poly_atom_z)+0.1)
fclose(file_id);

freez=ceil(poly_atom_number/2)
#create index
file_id=fopen("index.ndx","w");
fprintf(file_id,"[ System ]\n")
	for ii=1:poly_atom_number*mol_number
		fprintf(file_id,"%7d",ii)
		if(rem(ii,10)==0)
			fprintf(file_id,"\n");
		endif
	endfor
fprintf(file_id,"\n[ MF ]\n")
	for ii=1:mol_number
		fprintf(file_id,"%7d",freez+(ii-1)*poly_atom_number)
		if(rem(ii,10)==0)
			fprintf(file_id,"\n");
		endif
	endfor
	fprintf(file_id,"\n");
fclose(file_id);





my_max_list(L, M, I) :- nth0(I, L, M), \+ (member(E, L), abs(E) > abs(M)).

findPlace(Zs,RowSize,PosRow1,PosCol1) :- 
  my_max_list(Zs,Val,Pos),
  divmod(Pos,RowSize,PosRow0,PosCol0),
  PosRow1 is PosRow0+1,
  PosCol1 is PosCol0+1.

length_(Length, List) :- length(List, Length).

list2matrix(List, RowSize, Matrix) :-
    length(List, L),
    HowManyRows is L div RowSize,
    length(Matrix, HowManyRows),
    maplist(length_(RowSize), Matrix),
    append(Matrix, List).

extractCol(ColNumber, Matrix, Column) :-
    maplist(nth1(ColNumber), Matrix, Column).

extractSupports(ZsM,PosRow,PosCol,Support1,Support2) :-
    nth1(PosRow, ZsM, Support1),
    extractCol(PosCol, ZsM, Support2).

representList(L,RowSize,Support1,TransposedSupport2,ConstancyApp,ConstancyAddMeasurer,UnivarianceApp,UnivarianceAddMeasurer,A2) :-
    findPlace(L,RowSize,PosRow,PosCol),
    length(L,LLength),
    list2matrix(L,RowSize,Matrix),
    extractSupports(Matrix, PosRow, PosCol, Support1Raw, Support2Raw),
    vec_to_norm(Support1Raw, Support1Norm),
    vec_to_norm(Support2Raw, Support2Norm),
    ColSize is LLength div RowSize,
    list2matrix(Support1Raw,RowSize,Support1RawAsMatrix),
    list2matrix(Support2Raw,ColSize,Support2RawAsMatrix),
    matrix_div_scal(Support1RawAsMatrix, Support1Norm, Support1),
    matrix_div_scal(Support2RawAsMatrix, Support2Norm, Support2),
    transpose(Support1, TransposedSupport1),
    transpose(Support2, TransposedSupport2),
    matrix_multiply(Support2, Matrix, Temp),
    matrix_multiply(Temp, TransposedSupport1, ConstantComp),
    matrix_multiply(TransposedSupport2, Support1, Temp2),
    flatten(ConstantComp, ConstantCompAsList),
    max_list(ConstantCompAsList, ConstantCompAsElem),
    matrix_mult_scal(Temp2, ConstantCompAsElem, ConstancyApp),
    flatten(ConstancyApp, ConstancyAppAsList),
    vec_to_norm_sq(ConstancyAppAsList,ConstancyAppNorm),
    flatten(Matrix, MatrixAsList),
    vec_to_norm_sq(MatrixAsList,MatNorm),
    ConstancyAddMeasurer is ConstancyAppNorm / MatNorm,
    matrix_multiply(Matrix, TransposedSupport1, Temp3),
    matrix_mult_scal(TransposedSupport2, ConstantCompAsElem, Temp4),
    matrix_diff(Temp3, Temp4, Temp5),  % Temp5 is a(1)
    transpose(Matrix, TransposedMatrix),
    matrix_multiply(TransposedMatrix, TransposedSupport2, Temp6),
    matrix_mult_scal(TransposedSupport1, ConstantCompAsElem, Temp7),
    matrix_diff(Temp6, Temp7, Temp8),  % Temp8 is a(2),
    matrix_multiply(Temp5, Support1, Temp9),
    transpose(Temp8, Temp10),
    matrix_multiply(TransposedSupport2, Temp10, Temp11),
    matrix_sum(Temp11, Temp9, Temp12),  % a1vt+ua2T
    matrix_sum(ConstancyApp, Temp12, UnivarianceApp),
    flatten(UnivarianceApp, UnivarianceAppAsList),
    vec_to_norm_sq(UnivarianceAppAsList, UnivarianceAppNorm),
    UnivarianceAddMeasurer is UnivarianceAppNorm / MatNorm,
    matrix_diff(Matrix, UnivarianceApp, A2).
      
vec_to_norm_sq(L,Sum) :- foldl([X,FL,TR]>>(TR is X*X+FL),L,0,Sum).

vec_to_norm(L,Norm) :- 
     vec_to_norm_sq(L,Sum),
     Norm is sqrt(Sum).


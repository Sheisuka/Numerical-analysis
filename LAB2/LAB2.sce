// Точные методы решения СЛУ. Сравнение метода вращений и LU разложения

function[N, M, A, b] = read_data(filepath)
    f=mopen(filepath, 'r');
    N=mfscanf(f,'%d');
    M=mfscanf(f,'%d');
    for i=1:N
        for j=1:M
            A(i,j)=mfscanf(f,'%lg');
        end
    end
    for i=1:N
        b(i)=mfscanf(f,'%lg');
    end
    mclose(f);
endfunction

function[c, s] = get_cs(a1, a2)
     c = a1 / sqrt(a1 ^ 2 + a2 ^ 2);
     s = a2 / sqrt(a1 ^ 2 + a2 ^ 2);
endfunction
    
function[] = rotation_method(N, A)
    disp("Расшренная матрица");
    disp(A);
    for i = 1:N
        for j = i+1:N
            [c, s] = get_cs(A(i, i), A(j, i));
            Ai = A(i, :) .* c + A(j, :) .* s;
            Aj = A(j, :) .* c - A(i, :) .* s;
            A(i, :) = Ai;
            A(j, :) = Aj;
        end
    end
    disp("Привели к ступенчатому виду:");
    disp(A);
    
    x = zeros(1, N);
    for i = N:-1:1
        s = A(i, N + 1);
        for j = i+1:N
            s = s - A(i, j) * x(j);
        end;
        x(i) = s / A(i, i);
    end
    
    disp("Вычислили x");
    disp(x);
endfunction

[N, M, A, b] = read_data("data.txt");
matr = cat(2, A, b);
rotation_method(N, matr);




// Точные методы решения СЛУ. Сравнение метода вращений и LU разложения// Вариант 12

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
    
function[x] = rotation_method(N, A)
    for i = 1:N
        for j = i+1:N
            [c, s] = get_cs(A(i, i), A(j, i));
            Ai = A(i, :) .* c + A(j, :) .* s;
            Aj = A(j, :) .* c - A(i, :) .* s;
            A(i, :) = Ai;
            A(j, :) = Aj;
        end
    end
    
    x = zeros(1, N);
    for i = N:-1:1
        s = A(i, N + 1);
        for j = i+1:N
            s = s - A(i, j) * x(j);
        end;
        x(i) = s / A(i, i);
    end
endfunction

function[x] = LU_method(N, A, b)
    L = eye(N, N);
    U = zeros(N, N);
    
    for i = 1:N
        for j = i:N
            if i == 1
                s = 0;
            else
                s = L(i,1:i-1) * U(1:i-1,j);
            end
            U(i,j) = A(i,j) - s;
        end
        for j = i+1:N
            if i == 1
                s = 0;
            else
                s = L(j,1:i-1) * U(1:i-1,i);
            end
            L(j,i) = (A(j,i) - s) / U(i,i);
        end
    end
    
    y = zeros(1, N);
    for i = 1:N
        s = 0;
        for j = 1:i
            s = s + L(i, j) * y(j);
        end
        y(i) = b(i) - s;
    end
    
    x = zeros(1, N);
    for i = N:-1:1
        s = y(i);
        for j = i+1:N
            s = s - U(i, j) * x(j);
        end;
        x(i) = s / U(i, i);
    end
    
    disp(x);
endfunction


[N, M, A, b] = read_data("data.txt");
matr = cat(2, A, b);

disp("Метод вращений:")
x = rotation_method(N, matr);
disp("Найденное решение:", x);
r = A * x' - b;
disp("Вектор невязки:", r');
disp("Мера вектора невязки:", sqrt(r' * r));

disp("LU разложение:");
x = LU_method(N, A, b);
disp("Найденное решение:", x);
r = A * x' - b;
disp("Вектор невязки:", r');
disp("Мера вектора невязки:", sqrt(r' * r));

disp("Мера вектора невязки полученная при помощи решения с LU разложением получилась меньше")


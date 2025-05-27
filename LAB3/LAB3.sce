// Обусловленность СЛАУ
// Матрицы из прошлой лаб. работы - вариант 12

function[N, A, b] = read_data(filepath)
    f=mopen(filepath, 'r');
    N=mfscanf(f,'%d');
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
endfunction

function[A_h, b_h, x_h] = hilbert_matrix(N)
    A_h = zeros(N, N);
    x_h = ones(1, N);
    b_h = zeros(1, N);
    for i = 1:N
        for j = 1:N
            A_h(i, j) = 1 / (i + j - 1);
        end
        b_h(i) = A_h(i, :) * x_h';
    end
endfunction

[N] = read_data("data.txt"); // Считываем размерность, матрицу коэф., вектор св. членов 
[A_h, b_h, x_h] = hilbert_matrix(N); //  Генирируем матрицу Гильберта, решение и вектор св. членов
A_b = cat(2, A_h, b_h'); // Конкатенируем сгенерированные матрицу коэф. и вектор св. членов для решения методом вращений 
disp(A_b);
x_rot = rotation_method(N, A_b); // Ищем решение методов вращений
disp(x_rot);
x_lu = LU_method(N, A_h, b_h); // Ищем решение с помощью LU разложения
disp(x_lu);

r_rot = A_h * x_rot' - b_h'; // Вектор невязки для решения методом вращений
disp(r_rot);
r_rot_norm = sqrt(r_rot' * r_rot);
disp(r_rot_norm);
r_lu = A_h * x_lu' - b_h'; // Вектор невязки для решения с LU
disp(r_lu);
r_lu_norm = sqrt(r_lu' * r_lu);
disp(r_lu_norm);

cond_number_h = norm(A_h, "fro") * norm(inv(A_h), "fro"); // Число обусл. матрицы Гильберта
disp(cond_number_h);


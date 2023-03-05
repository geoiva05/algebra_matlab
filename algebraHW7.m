%% Задача 1 
% Решите систему методом сингулярного разложения:

A = [1 1 1; 1 3 1; 1 1 3];
b = [2; 4; 0];
[U, S, V] = svd(A);
x = U * (S ^ -1) * V' * b;
disp(x);

%% Задача 2 
% Решите систему из п. 1 методом разложения Холецкого. 

A = [1 1 1; 1 3 1; 1 1 3];
b = [2; 4; 0];
L = chol(A, 'lower');
y = L \ b;
x = L' \ y;
disp(x);

%% Задача 3 
% Напишите алгоритм итерационного метода Ричардсона 
% (см. источник 1, стр. 130, или слайд 4) и 
% решите с его помощью систему из пункта 1.

A = [1 1 1; 1 3 1; 1 1 3];
b = [2; 4; 0];
tau = 0.1;
x = [0; 0; 0];
n = 250;
for i=1:n
    r = b - A * x;
    x = x + r * tau;
end
disp(x)

%% Задача 4
% Напишите алгоритм метода простой итерации 
% (см. стр. 132 источника 1 или слайды 5-6) 
% и решите с его помощью систему из пункта 1.

A = [6 4 0; 1 3 1; 1 1 3];
b = [16; 4; 0];

x0 = [0; 0; 0];
n = 2000;
eps = 0.0001;

for i = 1:length(b)
    for j = 1:length(b)
        beta(i) = b(i) /A(i, i);
        if (i == j)
            newa(i, j) = 0;
        else 
            newa(i, j) = -A(i, j) / A(i, i);
        end
    end 
end

x1 = x0;
ncount = 0;
beta = beta';

while true
    ncount = ncount + 1;
    x1 = beta + newa * x0;
    max = abs(x0(1) - x1(1));
    for i=2:length(x0)
        if (abs(x0(i) - x1(i)) > max)
            max = abs(x0(i) - x1(i));
        end
    end
    if (max < eps) || (ncount > n)
        x = x1;
        disp(x)
        break;
    else 
        x0 = x1;
    end 
end 

%% Задача 5 
% Напишите алгоритм итерационного метода Гаусса-Зейделя 
% (см. источник 1, стр. 135 или слайды 7-8) и
% решите с его помощью систему из пункта 1.

A = [1 1 1; 1 3 1; 1 1 3];
b = [2; 4; 0];

x0 = [0; 0; 0];
n = 2000;
eps = 0.0001;
F = A' * A;
H = A' * b;

for i=1:length(b)
    for j=1:length(b)
        beta(i) = H(i) / F(i, i);
        if i == j
            newa(i, j) = 0;
        else 
            newa(i, j) = -F(i, j) / F(i, i);
        end
    end
end

x1 = x0;
ncount = 0;
beta = beta';

while true
    ncount  = ncount + 1;
    for i=1:length(b)
        s = 0;
        for j=1:length(b)
            s = s + newa(i, j) * x1(j);
        end
        x1(i) = beta(i) + s;
    end
    max = abs(x0(1) - x1(1));
    for i=2:length(x0)
        if abs(x0(i) - x1(i)) > max
            max = abs(x0(i) - x1(i));
        end
    end
    if (max < eps || ncount > n)
        x = x1;
        disp(x);
        break;
    else
        x0 = x1;
    end
end

%% Задача 6 Напишите алгоритм итерационного метода последовательной 
% верхней релаксации (SOR) (см. источник 1, стр. 136, или слайды 9-10) 
% и решите с его помощью систему из пункта 1.

A = [1 1 1; 1 3 1; 1 1 3];
b = [2; 4; 0];

x0 = [0; 0; 0];
n = 45;
eps = 0.00001;
F = A' * A;
H = A' * b;
w = 1.4;
for i=1:length(b)
    for j=1:length(b)
        beta(i) = H(i)/F(i, i);
        if i == j
            newa(i, j) = 0;
        else 
            newa(i, j) = -F(i, j) / F(i, i);
        end
    end
end

x1 = x0;
ncount = 0;
beta = beta';
while true 
    nctount = ncount + 1;
    for i=1:length(b)
        s = 0;
        for j=1:length(b)
            s = s + newa(i, j) * x1(j);
        end
        x1(i) = beta(i) + s + (w-1) * (beta(i) + s - x0(i));
    end
    max = abs(x0(1) - x1(1));
    for i=2:length(x0)
        if abs(x0(i) - x1(i)) > max
            max = abs(x0(i) - x1(i));
        end
    end
    if (max < eps || ncount > n)
        x = x1;
        disp(x);
        break;
    else
        x0 = x1;
    end
end

%% Задача 7 
% Напишите алгоритм итерационного метода сопряжённых градиентов
% (см. источник 1, стр. 181) и решите с его помощью систему из пункта 1.
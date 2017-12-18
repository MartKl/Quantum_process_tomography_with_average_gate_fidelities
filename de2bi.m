function y = de2bi(d,n)
% free and slow implementation f de2bi
y = zeros([1,n]);
for j = 1:n
    y(j) = rem(d,2);
    d = (d-y(j))./2;
end
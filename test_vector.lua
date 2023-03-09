require('vector')

l = vector.ListXf()
for i=1,10 do
	l:pushBack(i*i)
end
for i=0,9 do
	print(l[i])
end
for i=1,10 do
	local x = l:popBack()
	print(x)
end

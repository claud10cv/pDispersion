function extractQ(filename)
	Q = zeros(2, 0)
	f = open(filename)
	while !eof(f)
		line = readline(f)
		tok = split(line)
		if tok[1] == "Q"
			curr = 2
			nQ = floor(Int64, length(tok) / 2)
			for q in 1 : nQ
				app = [parse(Float64, tok[curr]), parse(Float64, tok[curr + 1])]
				curr += 2
				Q = hcat(Q, app)
			end
		end
	end
	close(f)
	Q
end

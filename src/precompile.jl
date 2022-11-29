@precompile_setup begin

    R = rand(10000)
    x = [1:100..., 1]
    y = [1:100..., NaN]
    Y = reshape(1:300,100,3)

    @precompile_all_calls begin

        histcounts(R, 0:0.1:1)
        histcounts(0.5:9.5,0.5:9.5, 0:10,0:10)
        nanbinmean(x, y, range(0,100,length=4))
        nanbinmean(1:100, Y, range(0,100,length=4))

        for T in (Float64,)
            for nd in 1:4
                A = ones(T, ntuple(i->10, nd))
                nansum(A)
                nanmean(A)
                nanstd(A)
                nanvar(A)
                nanminimum(A)
                nanmaximum(A)

                for d in 1:nd
                    nansum(A, dims=d)
                    nanmean(A, dims=d)
                    nanstd(A, dims=d)
                    nanvar(A, dims=d)
                    nanminimum(A, dims=d)
                    nanmaximum(A, dims=d)
                end
                if nd > 1
                    for i = 2:nd
                        for j = 1:i-1
                            nansum(A, dims=(j,i))
                            nanmean(A, dims=(j,i))
                            nanstd(A, dims=(j,i))
                            nanvar(A, dims=(j,i))
                            nanminimum(A, dims=(j,i))
                            nanmaximum(A, dims=(j,i))
                        end
                    end
                end

            end
        end
    end
end

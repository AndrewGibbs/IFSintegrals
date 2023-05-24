using NBInclude#, Test

# macro no_error(ex)
#     quote
#         try
#             $(esc(ex))
#             true
#         catch
#             false
#         end
#     end
# end

eg_files = readdir("../examples")
nb_files_inds = Int64[]
for (index,file) ∈ enumerate(eg_files)
    if file[(end-5):end]==".ipynb"
        push!(nb_files_inds,index)
    end
end
nb_files = eg_files[nb_files_inds]

# nb_files = filter((1:length(eg_files))[file[(end-5):end]==".ipynb" for file ∈ eg_files],eg_files)
# nb_files = ["Cantor set example.ipynb","BEM 3D example.ipynb","Quadrature example 1.ipynb","Quadrature example 2.ipynb"]

@testset begin
    for file ∈ nb_files
        println("../examples/"*file)
        @test@no_error(@nbinclude("../examples/"*file))
    end
end
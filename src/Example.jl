module Example

importall Base

hello(who::ASCIIString) = "Hello, $who"
helloworld() = println(hello("World"))

domath(x::Number) = (x + 6)

end

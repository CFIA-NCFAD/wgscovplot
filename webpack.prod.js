const path = require("path")
const TerserPlugin = require("terser-webpack-plugin");
module.exports = {
    mode: "production",
    entry: {
        wgscovplot:'./src/index.js',
    },
    optimization:{
        minimize: true,
        minimizer: [new TerserPlugin({
            extractComments: false,
            terserOptions: {
                format: {
                  comments: false,
                },
            }
        })]
    },
    output:{
        path: path.resolve(__dirname, './wgscovplot/tmpl/js'),
        filename: "[name].prod.bundle.js",
        library: "wgscovplot",
        libraryTarget: "umd",
        globalObject: "this"
    },
    module: {
        rules:[{
                test:/\.js$/,
                exclude:[/(node_modules)/,/(wgscovplot)/],
                use: "babel-loader"
            },
            {
                test: /\.css$/,
                use: ['style-loader', 'css-loader']
            },
            {
                test: require.resolve("jquery"),
                loader: "expose-loader",
                options: {
                    exposes: ["$", "jQuery"]
                }
            }]
    }
};
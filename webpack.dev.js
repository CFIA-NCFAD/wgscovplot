const path = require("path")
module.exports = {
    mode: "development",
    entry: {
        wgscovplot:'./src/index.js',
    },
    output:{
        path: path.resolve(__dirname, './wgscovplot/tmpl/js'),
        filename: "[name].dev.bundle.js",
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
        }],
    },
};
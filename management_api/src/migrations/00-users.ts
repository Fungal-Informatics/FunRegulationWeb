import { MigrationDefinitionWithName } from "umzug/lib/migrationsList";
import { zap } from "../helpers/sql";

const migration: MigrationDefinitionWithName = {
	name: "00-users",
	async up() {
		await zap.transaction.single/*sql*/`
            create table users (
                id uuid primary key,
                email text not null unique,
                login text not null, 
                password text not null,
                name text not null,
                "createdAt" timestamptz not null,
                "updatedAt" timestamptz not null,
				"isArchived" boolean not null
            );
        `;
	},
	async down() {
		await zap.transaction.single/*sql*/`
            drop table users;
        `;
	},
};

export default migration;
